%% Successive Convexification with Ellipsoidal Keep-Out Zone (KOZ)
% - Dynamics: Gim–Alfriend STM (J2) in LVLH (RTN) coordinates
% - Objective: minimize sum ||u_k||^2 + w_tr*sum ||x_k - x_bar_k||^2 + w_v*sum ||v_k||^2 + w_s*sum s_k
% - Constraints:
%     x_{k+1} = A_k x_k + B_k u_k + v_k
%     |u_i| <= u_max (component-wise)
%     Ellipsoidal KOZ avoidance (outside of ellipsoid) via linearization:
%       f(p) = (p-c)' Qinv (p-c) - 1 >= 0  (nonconvex)
%       Linearize at p_bar: f(p_bar) + grad' (p - p_bar) + s >= 0 (convex)
%     s_k >= 0 are slacks, penalized with w_s

clear; clc; close all;

%% Problem setup
% Physical parameters
params = struct();
params.Req = 6378.137e3;
params.mu = 3.986004415e14;
params.J2 = 1082.63e-6;
params.tol = 1e-12;
params.safetyAltitude = 100e3;

% Chief orbit
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

% Mission
mission = struct();
mission.u_max = 0.1; % [m/s^2] component-wise bound

% Time grid (length 60, interval 1 s)
dt = 1;  % [s]
N = 60;  % nodes
time_vec = linspace(0, dt*(N-1), N);

% Start/Goal override for SCP (meters): [x vx y vy z vz]
RelInitState  = [0; 0;  50; 0; 0; 0];
RelFinalState = [0; 0; -50; 0; 0; 0];
fprintf('Start (m): %s\n', sprintf('% .3f ', RelInitState));
fprintf('Goal  (m): %s\n', sprintf('% .3f ', RelFinalState));
fprintf('Init range = %.3f km, Final range = %.3f m\n', ...
    norm(RelInitState([1,3,5]))/1000, norm(RelFinalState([1,3,5])));

% Build STM
stm = local_computeSTM(params, chief, time_vec, RelInitState);

% KOZ definition (edit these)
koz = struct();
% Ellipsoid center (m) in LVLH [x;y;z]
koz.c = [0; 0; 0];
% Semi-axes (m) [ax, ay, az] from art-aeroconf24 (width/length/height mapping)
% height_ss+15 = 60, length_ss+20 = 94, width_ss+15 = 123
koz.axes = [60; 94; 123];
koz.use_ellipsoid = false; % disable the ellipsoidal KOZ (use implicit surfaces only)
% Apply KOZ only in a mid-horizon time window (fraction of [0,1])
koz.window = [0.0, 1.0]; % active across the whole horizon
% Guard margin: enforce f_lin >= margin (>0 pushes slightly outside)
koz.margin = 0.2;
% Optional spatial gate (activate KOZ only if |y|<=gate_T)
koz.gate_T = inf; % no gate; enforce everywhere

% Optional implicit-function KOZ list (union of surfaces)
% F_i(x,y,z) defined by external expressions: boundary when F_i-1=0.
% Outside should satisfy sign_i*(F_i-1) >= 0.
koz_list = {};
try
    expr_dir = 'C:\Users\98kim\Desktop\Acta-Astronautica\Funcs_ISS_expr';
    if exist(expr_dir,'dir')
        d = dir(fullfile(expr_dir,'*_*_expr.txt')); % use all available expressions
        take = numel(d);
        for i=1:take
            expr_path = fullfile(expr_dir, d(i).name);
            [Fh, Gh] = local_make_F_and_grad_from_expr(expr_path);
            % sign=+1 assumes F-1>=0 is outside. adjust if needed.
            koz_list{end+1} = struct('F', Fh, 'G', Gh, 'sign', +1, 'window', [0.4 0.85], 'margin', 0.02); %#ok<AGROW>
        end
    end
catch
end

% attach to koz structure with eval closures (h,grad)
koz.list = {};
for i=1:numel(koz_list)
    L = koz_list{i};
    if ~isempty(L.G)
        eval_handle = @(p) koz_implicit_eval_sym(L.F, L.G, p, L.sign);
    else
        eval_handle = @(p) koz_implicit_eval(L.F, p, L.sign);
    end
    koz.list{end+1} = struct('eval', eval_handle, 'window', L.window, 'margin', L.margin); %#ok<AGROW>
end
koz.Qinv = diag(1./(koz.axes.^2)); % shape matrix inverse

% SCP weights
weights = struct();
weights.R   = 1e4;   % control effort
weights.w_tr= 1e-5;  % trust region (state) penalty (looser movement)
weights.w_v = 1e2;   % virtual control penalty
weights.w_s = 1e8;   % KOZ slack penalty (much stronger)

max_iter   = 100;     % allow more iterations
min_iter   = 5;       % require at least this many iterations
tol_rel    = 5e-4;    % stricter state change tolerance
tol_cost   = 1e-4;    % relative cost improvement tolerance
tol_viol   = 1e-3;    % max KOZ violation tolerance (in f-units)
tol_vctrl  = 1e-3;    % virtual control norm tolerance
under_relax= 0.5;     % 0<beta<=1, step under-relaxation

% Initial guess: always solve an obstacle-free Q-only QP, then fallback to linear
Xbar = zeros(6,N); Ubar = zeros(3,N-1);
used_qp_seed = false;
thisdir = fileparts(mfilename('fullpath'));
try
    [XS,US] = local_qonly_seed(stm, N, RelInitState, RelFinalState, mission, dt);
    if all(isfinite(XS(:))) && all(isfinite(US(:))) && size(XS,2)==N && size(US,2)==N-1
        % Ensure exact boundary consistency
        XS(:,1)  = RelInitState;
        XS(:,end)= RelFinalState;
        Xbar = XS; Ubar = US; used_qp_seed = true;
        % Persist for reference with metadata (only if needed later)
        try
            X_opt = Xbar; %#ok<NASGU>
            U_opt = Ubar; %#ok<NASGU>
            save(fullfile(thisdir,'qp_qonly_last.mat'), 'X_opt','U_opt','time_vec','RelInitState','RelFinalState');
        catch
        end
    end
catch
end
if ~used_qp_seed
    for k=1:N
        a = (k-1)/(N-1);
        Xbar(:,k) = (1-a)*RelInitState + a*RelFinalState;
    end
    Ubar(:) = 0;
end

% Preserve the pre-SCP initial guess for plotting (do not overwrite during SCP)
Xbar_initial = Xbar; %#ok<NASGU>
Ubar_initial = Ubar; %#ok<NASGU>

best = struct('cost', inf, 'X', Xbar, 'U', Ubar);
prev_cost = inf; prev_viol = inf; prev_vnorm = inf;

Epos = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];

fprintf('\nRunning SCP (implicit KOZ surfaces, no ellipsoid)...\n');
for j = 1:max_iter
    [H,f,Aeq,beq,Ai,bi,lb,ub,idx] = build_qp(stm, N, RelInitState, RelFinalState, mission, ...
        weights, Xbar, koz, Epos);
    opts = optimoptions('quadprog','Display','off','MaxIterations',3000, ...
        'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8);
    [z, fval, exitflag] = quadprog(H,f,Ai,bi,Aeq,beq,lb,ub,[],opts);
    if exitflag<=0
        fprintf('  iter %d: solver failed (flag=%d)\n', j, exitflag);
        break;
    end
    [X,U,V,S] = unpack(z, idx, N);

    % Diagnostics
    if isfield(koz,'list') && ~isempty(koz.list)
        viol = koz_violation_multi(X, koz, Epos);
    else
        viol = koz_violation(X, koz, Epos);
    end
    cost = fval;
    vnorm = norm(V(:));
    dv   = computeTotalDeltaV(U, dt);
    rel_dX = norm(X(:)-Xbar(:)) / max(1, norm(Xbar(:)));
    rel_dJ = (prev_cost - cost) / max(1, abs(prev_cost));
    max_viol = max(viol);
    min_viol = min(viol);
    mean_viol = mean(viol);
    fprintf('  iter %02d: J=%.6e, dV=%.3f, viol[min/mean/max]=[%.2e/%.2e/%.2e], ||V||=%.3e, rel|dX|=%.2e, rel dJ=%.2e, w_s=%.2e\n', ...
        j, cost, dv, min_viol, mean_viol, max_viol, vnorm, rel_dX, rel_dJ, weights.w_s);

    % Update best
    if cost < best.cost
        best.cost = cost; best.X = X; best.U = U; best.V = V; best.S = S;
    end
    % Ramp slack penalty if violation persists
    if max_viol > tol_viol
        weights.w_s = min(weights.w_s*2, 1e9);
    end
    % Convergence check (require multiple criteria and min iterations)
    if j >= min_iter && rel_dX < tol_rel && abs(rel_dJ) < tol_cost && max_viol < tol_viol && vnorm < tol_vctrl
        fprintf('Converged: relX=%.2e, relJ=%.2e, viol=%.2e, v=%.2e\n', rel_dX, rel_dJ, max_viol, vnorm);
        Xbar = X; Ubar = U; break;
    end
    % Under-relaxed reference update
    Xbar = (1-under_relax)*Xbar + under_relax*X;
    Ubar = (1-under_relax)*Ubar + under_relax*U;
    prev_cost = cost; prev_viol = max_viol; prev_vnorm = vnorm;
end

X = best.X; U = best.U;

%% Console output: initial vs final trajectories
try
    fprintf('\n==== Initial trajectory (pre-SCP) ====\n');
    local_print_traj_table(Xbar_initial, time_vec);
catch
    % if Xbar_initial missing, skip without failing
end
fprintf('\n==== Final trajectory (optimized) ====\n');
local_print_traj_table(X, time_vec);

%% Plot
figure('Position',[60 60 1200 700]);
subplot(2,2,1);
hold on; grid on; box on;
% If a pre-rendered KOZ surface figure exists, load and copy it here
% and clear koz.list to avoid re-rendering heavy implicit surfaces.
try
    fig_file = 'C:\\Users\\98kim\\Desktop\\SCP\\Rendezvous\\ISS_data3_mesh120_opac05_dull.fig';
    if false && exist(fig_file,'file')
        ws = warning; warning('off','all');
        fh = openfig(fig_file, 'invisible');
        warning(ws);
        axsrc = findobj(fh,'Type','Axes');
        if ~isempty(axsrc)
            axsrc = axsrc(1);
            copyobj(allchild(axsrc), gca);
        end
        close(fh);
    end
catch
end
% Draw high‑res implicit surfaces over trajectory (if defined)
try
    if isempty(findobj(gca,'Type','Patch')) && isfield(koz,'list') && ~isempty(koz.list)
        lo = [-100; -100; -100]; hi = [100; 100; 100]; % meters domain box
        local_draw_implicit_surfaces(gca, koz, lo, hi, 100);
    end
catch
end
% Plot initial (pre‑SCP) and final trajectories
if exist('Xbar_initial','var')
    h_init = plot3(Xbar_initial(1,:)/1e3, Xbar_initial(3,:)/1e3, Xbar_initial(5,:)/1e3,'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
else
    h_init = gobjects(0);
end
h_final = plot3(X(1,:)/1e3, X(3,:)/1e3, X(5,:)/1e3,'b-','LineWidth',2);
h_start = plot3(X(1,1)/1e3, X(3,1)/1e3, X(5,1)/1e3,'go','MarkerFaceColor','g');
h_goal  = plot3(X(1,end)/1e3, X(3,end)/1e3, X(5,end)/1e3,'ro','MarkerFaceColor','r');
axis equal; view(45,30);
xlim([-0.1 0.1]); ylim([-0.1 0.1]); zlim([-0.1 0.1]); % +/-100 m (km units)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('3D Trajectory (pre/post) + implicit KOZ');
% Build legend dynamically depending on whether surfaces exist
h_surfaces = findobj(gca,'Type','Patch');
if ~isempty(h_surfaces)
    legend([h_surfaces(1) h_init h_final h_start h_goal], ...
           {'Implicit surfaces','Initial traj','Final traj','Start','Goal'}, ...
           'Location','northwest');
else
    legend([h_init h_final h_start h_goal], ...
           {'Initial traj','Final traj','Start','Goal'}, ...
           'Location','northwest');
end

subplot(2,2,2);
stairs(time_vec(1:end-1)/3600, vecnorm(U,2,1)*1e3,'k-','LineWidth',1.5); grid on;
xlabel('Time [h]'); ylabel('|u| [mm/s^2]'); title('Control magnitude');

subplot(2,2,3);
rel = sqrt(sum(X([1,3,5],:).^2,1));
plot(time_vec/3600, rel,'m-','LineWidth',2); grid on; xlabel('Time [h]'); ylabel('Range [m]');

subplot(2,2,4);
dv = [0 cumsum(vecnorm(U,2,1)*dt)];
plot(time_vec/3600, dv,'g-','LineWidth',2); grid on; xlabel('Time [h]'); ylabel('Cum. dV [m/s]');

print(gcf,'scp_koz_results.png','-dpng','-r150');

% Also overlay on the original FIG and save as a separate PNG
try
    if false && exist('fig_file','var') && exist(fig_file,'file')
        local_overlay_on_fig(fig_file, Xbar_initial, X);
    end
catch
end

% Persist latest SCP result for downstream comparison
try
    save(fullfile(fileparts(mfilename('fullpath')),'scp_last.mat'), ...
         'X','U','koz','time_vec');
catch
end

%% Helpers
function [RelInitState, RelFinalState] = load_custom_states()
    here = fileparts(mfilename('fullpath'));
    cands = {fullfile(pwd,'custom_states.csv'), fullfile(here,'custom_states.csv'), fullfile(here,'..','custom_states.csv')};
    csvp = '';
    for i=1:numel(cands), if exist(cands{i},'file'), csvp=cands{i}; break; end, end
    if isempty(csvp)
        error('custom_states.csv not found. Create it with two rows [x vx y vy z vz].');
    end
    Xio = readmatrix(csvp);
    RelInitState = Xio(1,:).';
    RelFinalState = Xio(2,:).';
end

function stm = local_computeSTM(params, chief, time_vec, RelInitState)
    B = zeros(6,3); B(2,1)=1; B(4,2)=1; B(6,3)=1;
    B_full = zeros(6,6); B_full(:,1:3)=B;
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

function [H,f,Aeq,beq,Ai,bi,lb,ub,idx] = build_qp(stm, N, x0, xf, mission, weights, Xbar, koz, Epos)
    nx=6; nu=3; nv=6; ns=1; % one slack per step
    nX = N*nx; nU=(N-1)*nu; nV=(N-1)*nv; nS=N*ns;
    nvar = nX+nU+nV+nS;

    % Indices
    idx.X = @(k) ((k-1)*nx+1):(k*nx);
    idx.U = @(k) (nX+(k-1)*nu+1):(nX+k*nu);
    idx.V = @(k) (nX+nU+(k-1)*nv+1):(nX+nU+k*nv);
    idx.S = @(k) (nX+nU+nV+(k-1)*ns+1):(nX+nU+nV+k*ns);

    % Objective H,f
    H = sparse(nvar,nvar); f = zeros(nvar,1);
    % Control effort
    for k=1:N-1
        iu = idx.U(k); H(iu,iu) = H(iu,iu) + 2*weights.R*eye(nu);
    end
    % Virtual control
    for k=1:N-1
        iv = idx.V(k); H(iv,iv) = H(iv,iv) + 2*weights.w_v*eye(nv);
    end
    % Trust region around Xbar
    for k=1:N
        ix=idx.X(k); H(ix,ix) = H(ix,ix) + 2*weights.w_tr*eye(nx);
        f(ix) = f(ix) - 2*weights.w_tr*Xbar(:,k);
    end
    % Slack L1
    for k=1:N
        is = idx.S(k); f(is) = f(is) + weights.w_s;
    end

    % Equality: dynamics + initial + final
    nEq = nx + nx*(N-1) + nx;
    Aeq = sparse(nEq, nvar); beq = zeros(nEq,1);
    row=0;
    % initial x
    Aeq(row+(1:nx), idx.X(1)) = eye(nx); beq(row+(1:nx)) = x0; row=row+nx;
    % dynamics
    for k=1:N-1
        rows = row+(1:nx);
        Aeq(rows, idx.X(k+1)) = eye(nx);
        Aeq(rows, idx.X(k))   = -stm.Ak(:,:,k);
        Aeq(rows, idx.U(k))   = -stm.Bk(:,1:3,k);
        Aeq(rows, idx.V(k))   = -eye(nx);
        row=row+nx;
    end
    % final
    Aeq(row+(1:nx), idx.X(N)) = eye(nx); beq(row+(1:nx)) = xf; row=row+nx;

    % Inequalities
    % u bounds and KOZ linearization + s>=0
    % Additional implicit KOZ surfaces count
    expr_ns = 0;
    persistent KOZ_LIST_CACHE
    if isempty(KOZ_LIST_CACHE), KOZ_LIST_CACHE = []; end
    if isfield(koz,'list') && ~isempty(koz.list)
        expr_ns = numel(koz.list);
    end
    % Pre-count rows
    % u upper/lower + optional (ellipsoid KOZ + s) + (expr KOZ + s for each surface)
    nI = 2*(N-1)*nu;
    if koz.use_ellipsoid
        nI = nI + 2*N;
    end
    nI = nI + 2*N*expr_ns;
    Ai = sparse(nI, nvar); bi = zeros(nI,1);
    r=0;
    % |u| <= u_max
    for k=1:N-1
        for j=1:nu
            e = zeros(1,nu); e(j)=1;
            r=r+1; Ai(r, idx.U(k)) =  e; bi(r) =  mission.u_max;
            r=r+1; Ai(r, idx.U(k)) = -e; bi(r) =  mission.u_max;
        end
    end
    % KOZ (ellipsoid) disabled unless koz.use_ellipsoid = true
    if koz.use_ellipsoid
        k_active_start = max(1, floor(koz.window(1)*(N-1))+1);
        k_active_end   = min(N,  ceil(koz.window(2)*(N-1))+1);
        for k=1:N
            pbar = Epos*Xbar(:,k);
            g = 2*koz.Qinv*(pbar - koz.c); % gradient (3x1)
            fbar = (pbar-koz.c)'*koz.Qinv*(pbar-koz.c) - 1; % >=0 outside
            a = (g)';
            % Guard margin: require f_lin >= margin
            b = a*pbar - fbar + koz.margin;
            % a*Epos*x_k + s_k >= b  ->  -(a*Epos)*x_k - s_k <= -b
            r=r+1;
            % spatial gate: only force when |T|<=gate_T
            gate_on = (abs(pbar(2)) <= koz.gate_T);
            if k>=k_active_start && k<=k_active_end && gate_on
                Ai(r, idx.X(k)) = -(a*Epos); Ai(r, idx.S(k)) = -1; bi(r) = -b;
            else
                % deactivated: 0 <= large number
                % leave Ai row as zeros; set bi large
                bi(r) = 1e9;
            end
            % s_k >= 0  ->  -s_k <= 0
            r=r+1; Ai(r, idx.S(k)) = -1; bi(r) = 0;
        end
    end

    % Additional KOZ from implicit surfaces list
    if isfield(koz,'list') && ~isempty(koz.list)
        for k=1:N
            pbar = Epos*Xbar(:,k);
            for i=1:numel(koz.list)
                L = koz.list{i};
                k0 = max(1, floor(L.window(1)*(N-1))+1);
                k1 = min(N,  ceil(L.window(2)*(N-1))+1);
                r=r+1;
                if k>=k0 && k<=k1 && (abs(pbar(2))<=koz.gate_T)
                    [hbar,gbar] = L.eval(pbar);
                    a = (gbar)';
                    b = a*pbar - hbar + L.margin;
                    Ai(r, idx.X(k)) = -(a*Epos);
                    Ai(r, idx.S(k)) = -1; % reuse same slack per-time (conservative)
                    bi(r) = -b;
                else
                    bi(r) = 1e9;
                end
                r=r+1; Ai(r, idx.S(k)) = -1; bi(r) = 0;
            end
        end
    end

    % Bounds (none on x,v)
    lb = -inf(nvar,1); ub = inf(nvar,1);
    % Enforce hard constraints by fixing slack variables s_k = 0
    for k=1:N
        lb(idx.S(k)) = 0; %#ok<AGROW>
        ub(idx.S(k)) = 0; %#ok<AGROW>
    end
end

function [X,U,V,S] = unpack(z, idx, N)
    nx=6; nu=3; nv=6; ns=1;
    X = zeros(nx,N); U=zeros(nu,N-1); V=zeros(nv,N-1); S=zeros(ns,N);
    for k=1:N,   X(:,k) = z(idx.X(k)); end
    for k=1:N-1, U(:,k) = z(idx.U(k)); V(:,k)=z(idx.V(k)); end
    for k=1:N,   S(:,k) = z(idx.S(k)); end
end

function v = koz_violation(X, koz, Epos)
    N=size(X,2); v=zeros(1,N);
    for k=1:N
        p=Epos*X(:,k); v(k)=(p-koz.c)'*koz.Qinv*(p-koz.c) - 1; % >=0 outside
    end
end

function total_dV = computeTotalDeltaV(U_opt, dt)
    total_dV = sum(vecnorm(U_opt,2,1))*dt;
end

function local_overlay_on_fig(fig_file, Xbar_initial, X)
    try
        ws = warning; warning('off','all');
        fh = openfig(fig_file, 'invisible');
        warning(ws);
    catch
        return;
    end
    try
        axs = findobj(fh,'Type','Axes');
        if isempty(axs), close(fh); return; end
        ax = axs(1);
        axes(ax); hold(ax,'on');
        try
            if exist('Xbar_initial','var') && ~isempty(Xbar_initial)
                plot3(ax, Xbar_initial(1,:)/1e3, Xbar_initial(3,:)/1e3, Xbar_initial(5,:)/1e3, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5, 'DisplayName','Initial traj');
            end
        catch
        end
        plot3(ax, X(1,:)/1e3, X(3,:)/1e3, X(5,:)/1e3, 'b-', 'LineWidth',2, 'DisplayName','Final traj');
        plot3(ax, X(1,1)/1e3, X(3,1)/1e3, X(5,1)/1e3, 'go', 'MarkerFaceColor','g', 'DisplayName','Start');
        plot3(ax, X(1,end)/1e3, X(3,end)/1e3, X(5,end)/1e3, 'ro', 'MarkerFaceColor','r', 'DisplayName','Goal');
        try, legend(ax,'show','Location','northwest'); catch, end
        outp = fullfile(fileparts(mfilename('fullpath')), 'scp_koz_overlay.png');
        try, print(fh, outp, '-dpng','-r150'); catch, end
        try, close(fh); catch, end
    catch
        try, close(fh); catch, end
    end
end

function local_print_traj_table(X, time_vec)
    % Print a trajectory as a table: k, t[s], [x vx y vy z vz]
    if nargin<2 || isempty(time_vec)
        time_vec = 0:(size(X,2)-1);
    end
    N = size(X,2);
    fprintf(' k     t[s]           x           vx           y           vy           z           vz\n');
    for k=1:N
        fprintf('%2d %9.3f %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n', ...
            k-1, time_vec(k), X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), X(6,k));
    end
end

% draw_ellipsoid helper removed (ellipsoidal KOZ not visualized)

function [X_opt,U_opt] = local_qonly_seed(stm, N, x0, xf, mission, dt)
    % Q-only QP: minimize sum_k x_k' Qx x_k subject to dynamics and bounds
    nx=6; nu=3; nX=N*nx; nU=(N-1)*nu; nvar=nX+nU;
    % State cost weights (penalize position strongly; small on velocity)
    Qx = diag([1, 1e-3, 1, 1e-3, 1, 1e-3]);
    H = sparse(nvar,nvar); f=zeros(nvar,1);
    for k=1:N
        ix=(k-1)*nx+(1:nx);
        H(ix,ix) = H(ix,ix) + 2*Qx;
    end
    % Equality: initial + dynamics + final
    Aeq = sparse(nx + nx*(N-1) + nx, nvar); beq = zeros(size(Aeq,1),1);
    row=0;
    Aeq(row+(1:nx), 1:nx) = eye(nx); beq(row+(1:nx))=x0; row=row+nx;
    for k=1:N-1
        rows=row+(1:nx);
        Aeq(rows, k*nx+(1:nx)) = eye(nx);
        Aeq(rows, (k-1)*nx+(1:nx)) = -stm.Ak(:,:,k);
        Aeq(rows, nX+(k-1)*nu+(1:nu)) = -stm.Bk(:,1:3,k);
        row=row+nx;
    end
    Aeq(row+(1:nx), (N-1)*nx+(1:nx)) = eye(nx); beq(row+(1:nx))=xf;
    % Inequality: |u| <= u_max
    nI = 2*(N-1)*nu; Ai = sparse(nI,nvar); bi=zeros(nI,1); r=0;
    for k=1:N-1
        for j=1:nu
            e=zeros(1,nu); e(j)=1; r=r+1; Ai(r, nX+(k-1)*nu+(1:nu)) =  e; bi(r)= mission.u_max;
            r=r+1; Ai(r, nX+(k-1)*nu+(1:nu)) = -e; bi(r)= mission.u_max;
        end
    end
    % Solve
    opts = optimoptions('quadprog','Display','off','MaxIterations',3000,'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);
    z = quadprog(H,f,Ai,bi,Aeq,beq,[],[],[],opts);
    X_opt = reshape(z(1:nX), nx, N);
    U_opt = reshape(z(nX+1:end), nu, N-1);
end

function Fh = local_make_F_from_expr(expr_path)
    % Build a MATLAB anonymous function @(x,y,z) from a Python-style expr file
    txt = strtrim(fileread(expr_path));
    % Replace Python power operator ** with MATLAB .^ (scalars fine)
    txt = regexprep(txt, '\*\*', '.^');
    % Build function string
    fstr = sprintf('@(x,y,z) (%s)', txt);
    Fh = str2func(fstr);
end

function [Fh, Gh] = local_make_F_and_grad_from_expr(expr_path)
    % Try to build F and exact gradient using Symbolic Math; fallback to F only
    Fh = []; Gh = [];
    try
        txt = strtrim(fileread(expr_path));
        txt = regexprep(txt, '\*\*', '.^'); % Python power to MATLAB
        % symbolic variables
        syms xs ys zs real
        fsym = str2sym(strrep(strrep(strrep(txt,'x','xs'),'y','ys'),'z','zs'));
        Fh  = matlabFunction(fsym, 'Vars', [xs, ys, zs]);
        dfx = matlabFunction(diff(fsym, xs), 'Vars', [xs, ys, zs]);
        dfy = matlabFunction(diff(fsym, ys), 'Vars', [xs, ys, zs]);
        dfz = matlabFunction(diff(fsym, zs), 'Vars', [xs, ys, zs]);
        Gh = struct('dx', dfx, 'dy', dfy, 'dz', dfz);
    catch
        % Fallback to anonymous-only (numeric gradient will be used)
        try
            Fh = local_make_F_from_expr(expr_path);
        catch
            Fh = [];
        end
        Gh = [];
    end
end

function [h, grad] = koz_implicit_eval(F, p, sgn)
    % Evaluate h = sgn*(F(p)-1) and its gradient via finite differences
    x=p(1); y=p(2); z=p(3);
    f0 = F(x,y,z);
    epsv = 1e-4;
    fxp = F(x+epsv,y,z); fxm = F(x-epsv,y,z);
    fyp = F(x,y+epsv,z); fym = F(x,y-epsv,z);
    fzp = F(x,y,z+epsv); fzm = F(x,y,z-epsv);
    dfdx = (fxp - fxm)/(2*epsv);
    dfdy = (fyp - fym)/(2*epsv);
    dfdz = (fzp - fzm)/(2*epsv);
    h = sgn*(f0 - 1);
    grad = sgn*[dfdx; dfdy; dfdz];
end

function [h, grad] = koz_implicit_eval_sym(F, G, p, sgn)
    % Evaluate h and exact gradient using provided gradient functions G
    x=p(1); y=p(2); z=p(3);
    f0 = F(x,y,z);
    dfdx = G.dx(x,y,z);
    dfdy = G.dy(x,y,z);
    dfdz = G.dz(x,y,z);
    h = sgn*(f0 - 1);
    grad = sgn*[dfdx; dfdy; dfdz];
end

function local_draw_implicit_surfaces(ax, koz, lo, hi, res)
    axes(ax); %#ok<LAXES>
    colors = lines(max(1, numel(koz.list)));
    xv = linspace(lo(1),hi(1),res);
    yv = linspace(lo(2),hi(2),res);
    zv = linspace(lo(3),hi(3),res);
    [Xg,Yg,Zg] = meshgrid(xv,yv,zv);
    for i=1:numel(koz.list)
        L = koz.list{i};
        V = arrayfun(@(x,y,z) local_eval_h_val(L, x, y, z), Xg, Yg, Zg);
        vmin = min(V(:)); vmax = max(V(:));
        if vmin<=0 && vmax>=0
            try
                p = patch(isosurface(Xg/1e3,Yg/1e3,Zg/1e3,V,0));
                c = colors(1+mod(i-1,size(colors,1)),:);
                set(p,'FaceColor',c,'EdgeColor','none','FaceAlpha',0.25);
                isonormals(Xg/1e3,Yg/1e3,Zg/1e3,V,p);
            catch
            end
        end
    end
end

function h = local_eval_h_val(L, x, y, z)
    try
        [hval,~] = L.eval([x;y;z]);
        h = hval;
    catch
        [hval,~] = L.eval([x;y;z]);
        h = hval;
    end
end
function viol = koz_violation_multi(X, koz, Epos)
    N=size(X,2); viol=zeros(1,N);
    for k=1:N
        p = Epos*X(:,k);
        if isfield(koz,'list') && ~isempty(koz.list)
            vals = zeros(1,numel(koz.list));
            for i=1:numel(koz.list)
                [hval,~] = koz.list{i}.eval(p);
                vals(i)=hval;
            end
            viol(k) = min(vals);
        else
            viol(k) = 0;
        end
    end
end
