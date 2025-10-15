% Run QCQP-like rendezvous using custom start/end from CSV.
% If Gurobi is unavailable, falls back to fmincon with the same
% quadratic objective and second-order norm constraints.

clear; clc; close all;

% Inputs
csv_path = 'custom_states.csv';
N_override = []; % set numeric to force N; empty uses suggestion or default
dt = 60;         % [s] time step

% Read custom states: two rows [x vx y vy z vz] in meters
Xio = readmatrix(csv_path);
if size(Xio,1) < 2 || size(Xio,2) < 6
    error('custom_states.csv must have 2 rows and 6 columns');
end
RelInitState = Xio(1,:).';
RelFinalState = Xio(2,:).';

% Problem parameters (same as QCQP script defaults)
params = struct();
params.Req = 6378.137e3;
params.mu = 3.986004415e14;
params.J2 = 1082.63e-6;
params.tol = 1e-12;
params.safetyAltitude = 100e3;

chief = struct();
chief.altitude = 700e3;
chief.a = params.Req + chief.altitude;
chief.e = 0.001;
chief.i = deg2rad(98.2);
chief.RAAN = deg2rad(30);
chief.w = deg2rad(0);
chief.M = deg2rad(60);
chief.elements = [chief.a, chief.e, chief.i, chief.RAAN, chief.w, chief.M];
chief.n = sqrt(params.mu/chief.a^3);
chief.T = 2*pi/chief.n;

mission = struct();
% u_max pick conservative; adjust if needed
mission.u_max = 0.05; % [m/s^2]

% Build time grid
if isempty(N_override)
    % If a helper file with suggestion exists, we keep default N=101
    N = 101;
else
    N = N_override;
end
time_vec = linspace(0, dt*(N-1), N);

% Build STM (Ak,Bk)
stm = local_computeSTM(params, chief, time_vec, RelInitState);

% Solve via QCQP if Gurobi exists, otherwise fmincon fallback
if exist('gurobi','file')
    time_opts = struct('dt', dt);
    [sol, info] = solveQCQP_L2(stm, N, RelInitState, RelFinalState, mission, time_opts);
else
    fprintf('Gurobi not available. Using fmincon fallback.\n');
    [sol, info] = solveQCQP_L2_fmincon(stm, N, RelInitState, RelFinalState, mission, dt);
end

fprintf('\nStatus: %s\nSolve time: %.3f s\nEnergy: %.6f\nTotal dV: %.3f m/s\nFinal pos err: %.2e m\nFinal vel err: %.2e m/s\n',...
    info.status, info.solve_time, info.energy, info.total_dV, info.final_pos_error, info.final_vel_error);

% Save figures
local_visualize_basic(sol, time_vec, norm(RelFinalState([1,3,5])));
saveas(gcf, 'qcqp_custom_results.png');

% ---- Helper: fmincon fallback ----
function [sol, info] = solveQCQP_L2_fmincon(stm, N, RelInitState, RelFinalState, mission, dt)
    n_states = 6; n_controls = 3;
    n_x = N*n_states; n_u = (N-1)*n_controls; n_vars = n_x + n_u;

    % Equality constraints Aeq*z = beq
    Aeq = sparse(n_states + n_states*(N-1) + n_states, n_vars);
    beq = zeros(size(Aeq,1),1);
    row = 0;
    % initial
    Aeq(1:n_states,1:n_states) = eye(n_states); beq(1:n_states) = RelInitState; row = n_states;
    % dynamics
    for k=1:N-1
        rows = row+(1:n_states);
        Aeq(rows, k*n_states+(1:n_states)) = eye(n_states);
        Aeq(rows, (k-1)*n_states+(1:n_states)) = -stm.Ak(:,:,k);
        Aeq(rows, n_x+(k-1)*n_controls+(1:n_controls)) = -stm.Bk(:,1:3,k);
        row = row + n_states;
    end
    % final
    rows = row+(1:n_states);
    Aeq(rows, (N-1)*n_states+(1:n_states)) = eye(n_states);
    beq(rows) = RelFinalState;

    % Objective: sum ||u_k||^2
    function J = obj(z)
        U = reshape(z(n_x+1:end), n_controls, N-1);
        J = sum(sum(U.^2));
    end

    % Nonlinear constraints: ||u_k||_2 <= u_max
    function [c,ceq] = nonlcon(z)
        U = reshape(z(n_x+1:end), n_controls, N-1);
        c = vecnorm(U,2,1) - mission.u_max; % each <= 0
        ceq = [];
    end

    % Initial guess: linear interp for X, zero U
    X0 = zeros(n_states,N);
    alpha = linspace(0,1,N);
    for k=1:N
        xk = (1-alpha(k))*RelInitState + alpha(k)*RelFinalState;
        X0(:,k) = xk;
    end
    U0 = zeros(n_controls,N-1);
    z0 = [X0(:); U0(:)];

    opts = optimoptions('fmincon','Algorithm','interior-point','Display','iter-detailed',...
        'SpecifyObjectiveGradient',false,'MaxIterations',300,'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-6);

    tic;
    [z_opt, fval, exitflag, output] = fmincon(@obj, z0, [],[], Aeq, beq, [],[], @nonlcon, opts);
    solve_time = toc;

    X_opt = reshape(z_opt(1:n_x), n_states, N);
    U_opt = reshape(z_opt(n_x+1:end), n_controls, N-1);

    energy = fval;
    total_dV = sum(vecnorm(U_opt,2,1))*dt;
    final_pos_error = norm(X_opt([1,3,5],end) - RelFinalState([1,3,5]));
    final_vel_error = norm(X_opt([2,4,6],end) - RelFinalState([2,4,6]));

    sol = struct('X', X_opt, 'U', U_opt);
    info = struct('energy', energy, 'total_dV', total_dV, 'final_pos_error', final_pos_error, ...
        'final_vel_error', final_vel_error, 'solve_time', solve_time, 'status', string(exitflag));
end

% ---- Helper: compute STM ----
function stm = local_computeSTM(params, chief, time_vec, RelInitState)
    B = zeros(6,3);
    B(2,1) = 1; B(4,2) = 1; B(6,3) = 1;
    B_full = zeros(6,6);
    B_full(:,1:3) = B;
    initStruct.params = {params.Req, params.mu, params.J2, params.tol, params.safetyAltitude};
    initStruct.maneuverParams = {10, B_full};
    initStruct.timeParams.t0 = time_vec(1);
    initStruct.timeParams.dt = mean(diff(time_vec));
    initStruct.timeParams.tf = time_vec(end);
    initStruct.initChiefDescription = 'Classical';
    initStruct.initDeputyDescription = 'Cartesian';
    initStruct.Elements = chief.elements;
    initStruct.RelInitState = RelInitState;
    initStruct.RelFinalState = [];
    stm = GimAlfriendSTM(initStruct);
    stm.makeDiscreteMatrices();
end

% ---- Helper: basic LVLH visualization ----
function local_visualize_basic(sol, time_vec, rf)
    X = sol.X; U = sol.U; dt = mean(diff(time_vec));
    figure('Position',[50 50 1200 800]);
    subplot(2,2,1);
    plot3(X(1,:)/1e3, X(3,:)/1e3, X(5,:)/1e3,'b-','LineWidth',2); grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]'); title('LVLH Trajectory'); view(45,30);
    subplot(2,2,2);
    rel_dist = sqrt(sum(X([1,3,5],:).^2,1));
    plot(time_vec/3600, rel_dist,'k-','LineWidth',2); hold on; yline(rf,'r--'); grid on;
    xlabel('Time [h]'); ylabel('Relative distance [m]'); title('Range over time');
    subplot(2,2,3);
    stairs(time_vec(1:end-1)/3600, vecnorm(U,2,1)*1e3,'m-','LineWidth',1.5); grid on;
    xlabel('Time [h]'); ylabel('|u| [mm/s^2]'); title('Control magnitude');
    subplot(2,2,4);
    dv_cum = [0, cumsum(vecnorm(U,2,1)*dt)]; % length N
    plot(time_vec/3600, dv_cum,'g-','LineWidth',2); grid on;
    xlabel('Time [h]'); ylabel('Cum. dV [m/s]'); title('Cumulative dV');
end
