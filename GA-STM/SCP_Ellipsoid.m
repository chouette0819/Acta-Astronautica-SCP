%% SCP with Ellipsoidal KOZ (Pipeline 1)
% - Read start/goal from CSV
% - Build QP initial guess (Q-only state cost)
% - Run SCP with ellipsoidal KOZ only (no implicit surfaces)

clear; clc; close all;

% Inputs
start_goal_csv = 'custom_states.csv'; % two rows: [x vx y vy z vz] (meters)

% Time grid
dt = 60;          % [s]
N = 101;          % nodes
time_vec = linspace(0, dt*(N-1), N);

% Physical parameters
params = struct();
params.Req = 6378.137e3; params.mu = 3.986004415e14; params.J2 = 1082.63e-6;
params.tol = 1e-12; params.safetyAltitude = 100e3;

% Chief orbit (dummy, only for STM)
chief = struct();
chief.altitude = 700e3; chief.a = params.Req + chief.altitude;
chief.e = 0.001; chief.i = deg2rad(98.2); chief.RAAN = deg2rad(30);
chief.w = deg2rad(0); chief.M = deg2rad(60);
chief.elements = [chief.a, chief.e, chief.i, chief.RAAN, chief.w, chief.M];

% Load start/goal
[RelInitState, RelFinalState] = load_custom_states_fallback(start_goal_csv);
fprintf('Start (m): %s\n', sprintf('% .3f ', RelInitState));
fprintf('Goal  (m): %s\n', sprintf('% .3f ', RelFinalState));

% Build STM
stm = local_computeSTM(params, chief, time_vec, RelInitState);

% Ellipsoidal KOZ (edit as needed)
koz = struct();
koz.c = [0; 0; 0];            % center [m]
koz.axes = [60; 94; 123];     % semi-axes [m]
koz.margin = 0.0;             % guard margin on linearization
koz.window = [0.0, 1.0];      % active across whole horizon
koz.gate_T = inf;             % no T-gate
koz.Qinv = diag(1./(koz.axes.^2));

% Weights
weights = struct('R',1e4,'w_tr',1e-3,'w_v',1e2,'w_s',1e6);

% Solve
max_iter=50; min_iter=3; tol_rel=1e-3; tol_cost=1e-4; tol_viol=1e-3; tol_vctrl=1e-3; beta=0.5;

% Initial guess via Q-only QP
[Xbar,Ubar] = local_qonly_seed(stm, N, RelInitState, RelFinalState, struct('u_max',0.1), dt);

% SCP loop
Epos = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
best = struct('cost', inf, 'X', Xbar, 'U', Ubar);
prev_cost = inf;
fprintf('\nRunning SCP (Ellipsoid KOZ) ...\n');
for j=1:max_iter
    [H,f,Aeq,beq,Ai,bi,lb,ub,idx] = build_qp_ellipsoid(stm,N,RelInitState,RelFinalState,struct('u_max',0.1),weights,Xbar,koz,Epos);
    opts = optimoptions('quadprog','Display','off','MaxIterations',3000,'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8);
    [z, fval, exitflag] = quadprog(H,f,Ai,bi,Aeq,beq,lb,ub,[],opts);
    if exitflag<=0, fprintf('  iter %02d: solver failed (flag=%d)\n', j, exitflag); break; end
    [X,U,V,S] = unpack(z, idx, N);
    % diag
    vnorm = norm(V(:)); dv = sum(vecnorm(U,2,1))*dt;
    viol = koz_violation_ellipsoid(X, koz, Epos);
    rel_dX = norm(X(:)-Xbar(:))/max(1,norm(Xbar(:)));
    rel_dJ = (prev_cost - fval)/max(1,abs(prev_cost));
    fprintf('  iter %02d: J=%.6e, dV=%.3f, min/max viol=[%.2e/%.2e], ||V||=%.2e, rel|dX|=%.2e, rel dJ=%.2e\n', ...
        j, fval, dv, min(viol), max(viol), vnorm, rel_dX, rel_dJ);
    if fval<best.cost, best.cost=fval; best.X=X; best.U=U; best.V=V; best.S=S; end
    if j>=min_iter && rel_dX<tol_rel && abs(rel_dJ)<tol_cost && max(viol)<tol_viol && vnorm<tol_vctrl
        fprintf('Converged at iter %d\n', j); Xbar=X; Ubar=U; break; end
    Xbar = (1-beta)*Xbar + beta*X; Ubar = (1-beta)*Ubar + beta*U; prev_cost=fval;
end

X=best.X; U=best.U;

% Plot
figure('Position',[60 60 1100 700]); subplot(2,2,1);
plot3(X(1,:)/1e3,X(3,:)/1e3,X(5,:)/1e3,'b-','LineWidth',2); hold on; grid on;
[xe,ye,ze] = ellipsoid(koz.c(1)/1e3,koz.c(2)/1e3,koz.c(3)/1e3,koz.axes(1)/1e3,koz.axes(2)/1e3,koz.axes(3)/1e3,20);
surf(xe,ye,ze,'FaceAlpha',0.15,'EdgeColor','none','FaceColor',[1 0 0]);
axis equal; view(45,30);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]'); title('SCP (Ellipsoid KOZ)');

% Save
try
    save(fullfile(fileparts(mfilename('fullpath')),'scp_ellipsoid_last.mat'),'X','U','time_vec');
catch, end

%% Helpers (subset)
function [RelInitState, RelFinalState] = load_custom_states_fallback(csvp)
    if exist(csvp,'file')
        Xio = readmatrix(csvp); RelInitState = Xio(1,:).'; RelFinalState = Xio(2,:).';
    else
        error('Start/goal file not found: %s', csvp);
    end
end

function stm = local_computeSTM(params, chief, time_vec, RelInitState)
    B = zeros(6,3); B(2,1)=1; B(4,2)=1; B(6,3)=1; B_full=zeros(6,6); B_full(:,1:3)=B;
    initStruct.params = {params.Req, params.mu, params.J2, params.tol, params.safetyAltitude};
    initStruct.maneuverParams = {10, B_full};
    initStruct.timeParams.t0 = time_vec(1); initStruct.timeParams.dt = mean(diff(time_vec)); initStruct.timeParams.tf = time_vec(end);
    initStruct.initChiefDescription='Classical'; initStruct.initDeputyDescription='Cartesian'; initStruct.Elements = chief.elements; initStruct.RelInitState = RelInitState; initStruct.RelFinalState = [];
    stm = GimAlfriendSTM(initStruct); stm.makeDiscreteMatrices();
end

function [H,f,Aeq,beq,Ai,bi,lb,ub,idx] = build_qp_ellipsoid(stm,N,x0,xf,mission,weights,Xbar,koz,Epos)
    nx=6; nu=3; nv=6; ns=1; nX=N*nx; nU=(N-1)*nu; nV=(N-1)*nv; nS=N*ns; nvar=nX+nU+nV+nS;
    idx.X=@(k)((k-1)*nx+1):(k*nx); idx.U=@(k)(nX+(k-1)*nu+1):(nX+k*nu); idx.V=@(k)(nX+nU+(k-1)*nv+1):(nX+nU+k*nv); idx.S=@(k)(nX+nU+nV+(k-1)*ns+1):(nX+nU+nV+k*ns);
    H=sparse(nvar,nvar); f=zeros(nvar,1);
    for k=1:N-1, H(idx.U(k),idx.U(k))=H(idx.U(k),idx.U(k))+2*weights.R*eye(nu); end
    for k=1:N-1, H(idx.V(k),idx.V(k))=H(idx.V(k),idx.V(k))+2*weights.w_v*eye(nv); end
    for k=1:N, H(idx.X(k),idx.X(k))=H(idx.X(k),idx.X(k))+2*weights.w_tr*eye(nx); f(idx.X(k))=f(idx.X(k))-2*weights.w_tr*Xbar(:,k); end
    for k=1:N, f(idx.S(k))=f(idx.S(k))+weights.w_s; end
    Aeq=sparse( N*nx + nx , nvar); beq=zeros(size(Aeq,1),1); row=0;
    Aeq(1:nx, idx.X(1))=eye(nx); beq(1:nx)=x0; row=row+nx;
    for k=1:N-1, rows=row+(1:nx); Aeq(rows,idx.X(k+1))=eye(nx); Aeq(rows,idx.X(k))=-stm.Ak(:,:,k); Aeq(rows,idx.U(k))=-stm.Bk(:,1:3,k); Aeq(rows,idx.V(k))=-eye(nx); row=row+nx; end
    Aeq(row+(1:nx), idx.X(N))=eye(nx); beq(row+(1:nx))=xf;
    % Inequalities: u box + ellipsoid + s>=0
    nI=2*(N-1)*nu + 2*N; Ai=sparse(nI,nvar); bi=zeros(nI,1); r=0;
    for k=1:N-1, for j=1:nu, e=zeros(1,nu); e(j)=1; r=r+1; Ai(r,idx.U(k))= e; bi(r)=mission.u_max; r=r+1; Ai(r,idx.U(k))=-e; bi(r)=mission.u_max; end, end
    k0=max(1,floor(koz.window(1)*(N-1))+1); k1=min(N,ceil(koz.window(2)*(N-1))+1);
    for k=1:N
        pbar=Epos*Xbar(:,k); g=2*koz.Qinv*(pbar-koz.c); fbar=(pbar-koz.c)'*koz.Qinv*(pbar-koz.c)-1; a=g'; b=a*pbar - fbar + koz.margin;
        r=r+1; if k>=k0 && k<=k1 && (abs(pbar(2))<=koz.gate_T), Ai(r,idx.X(k))=-(a*Epos); Ai(r,idx.S(k))=-1; bi(r)=-b; else, bi(r)=1e9; end
        r=r+1; Ai(r,idx.S(k))=-1; bi(r)=0;
    end
    lb=-inf(nvar,1); ub=inf(nvar,1);
end

function [X,U,V,S] = unpack(z,idx,N)
    nx=6; nu=3; nv=6; ns=1; nX=N*nx; nU=(N-1)*nu; nV=(N-1)*nv; nS=N*ns;
    X=zeros(nx,N); U=zeros(nu,N-1); V=zeros(nv,N-1); S=zeros(ns,N);
    for k=1:N, X(:,k)=z(idx.X(k)); end
    for k=1:N-1, U(:,k)=z(idx.U(k)); V(:,k)=z(idx.V(k)); end
    for k=1:N, S(:,k)=z(idx.S(k)); end
end

function viol=koz_violation_ellipsoid(X,koz,Epos)
    N=size(X,2); viol=zeros(1,N); for k=1:N, p=Epos*X(:,k); viol(k)=(p-koz.c)'*koz.Qinv*(p-koz.c)-1; end
end

function [X_opt,U_opt] = local_qonly_seed(stm, N, x0, xf, mission, dt)
    nx=6; nu=3; nX=N*nx; nU=(N-1)*nu; nvar=nX+nU; Qx=diag([1,1e-3,1,1e-3,1,1e-3]);
    H=sparse(nvar,nvar); f=zeros(nvar,1); for k=1:N, ix=(k-1)*nx+(1:nx); H(ix,ix)=H(ix,ix)+2*Qx; end
    Aeq=sparse(nx+nx*(N-1)+nx, nvar); beq=zeros(size(Aeq,1),1); row=0;
    Aeq(1:nx,1:nx)=eye(nx); beq(1:nx)=x0; row=row+nx;
    for k=1:N-1, rows=row+(1:nx); Aeq(rows,k*nx+(1:nx))=eye(nx); Aeq(rows,(k-1)*nx+(1:nx))=-stm.Ak(:,:,k); Aeq(rows,nX+(k-1)*nu+(1:nu))=-stm.Bk(:,1:3,k); row=row+nx; end
    Aeq(row+(1:nx),(N-1)*nx+(1:nx))=eye(nx); beq(row+(1:nx))=xf;
    nI=2*(N-1)*nu; Ai=sparse(nI,nvar); bi=zeros(nI,1); r=0; for k=1:N-1, for j=1:nu, e=zeros(1,nu); e(j)=1; r=r+1; Ai(r,nX+(k-1)*nu+(1:nu))= e; bi(r)=mission.u_max; r=r+1; Ai(r,nX+(k-1)*nu+(1:nu))=-e; bi(r)=mission.u_max; end, end
    opts=optimoptions('quadprog','Display','off','MaxIterations',3000); z=quadprog(H,f,Ai,bi,Aeq,beq,[],[],[],opts);
    X_opt=reshape(z(1:nX),nx,N); U_opt=reshape(z(nX+1:end),nu,N-1);
end

