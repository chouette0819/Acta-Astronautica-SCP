%% MIQP-based KOZ avoidance (L2 objective, same domain as SCP)
% - Dynamics: Gim–Alfriend STM (J2) in LVLH (RTN) coordinates
% - Objective (quadratic): minimize sum ||u_k||^2 + w_tr*sum ||x_k - x_bar_k||^2 + w_v*sum ||v_k||^2
% - Constraints:
%     x_{k+1} = A_k x_k + B_k u_k + v_k
%     |u_i| <= u_max (component-wise)
%     AABB KOZ (per implicit F file): outside-of-box via big-M disjunction with binaries
%
% Notes
% - Requires Gurobi MATLAB API on path (MIQP). If not found, errors out.
% - Time grid and start/goal are set to match SCP_KOZ_QP defaults.

clear; clc; close all;
% Solver selection: 'quadprog' (QP relaxation) or 'gurobi' (exact MIQP if available)
solver = 'quadprog';  % change to 'gurobi' to use Gurobi MIQP

%% Problem setup (match SCP_KOZ_QP)
params = struct();
params.Req = 6378.137e3; params.mu = 3.986004415e14; params.J2 = 1082.63e-6; params.tol = 1e-12; params.safetyAltitude = 100e3;

chief = struct();
chief.altitude = 700e3; chief.a = params.Req + chief.altitude;
chief.e = 0.001; chief.i = deg2rad(98.2); chief.RAAN = deg2rad(30);
chief.w = deg2rad(0); chief.M = deg2rad(60);
chief.elements = [chief.a, chief.e, chief.i, chief.RAAN, chief.w, chief.M];

mission = struct();
mission.u_max = 0.1; % [m/s^2] component-wise bound

% Time grid (same as SCP)
dt = 1;  % [s]
N = 60;  % nodes
time_vec = linspace(0, dt*(N-1), N);

% Start/Goal (same as SCP)
RelInitState  = [0; 0;  50; 0; 0; 0];
RelFinalState = [0; 0; -50; 0; 0; 0];

fprintf('Start (m): %s\n', sprintf('% .3f ', RelInitState));
fprintf('Goal  (m): %s\n', sprintf('% .3f ', RelFinalState));

% Build STM
stm = local_computeSTM(params, chief, time_vec, RelInitState);

% KOZ AABB options
koz = struct();
koz.expr_dir   = 'C:\\Users\\98kim\\Desktop\\Acta-Astronautica\\Funcs_ISS_expr';
koz.domain_lo  = [-100; -100; -100];  % [m]
koz.domain_hi  = [ 100;  100;  100];  % [m]
koz.grid_res   = 60;
koz.max_expr   = 12;
koz.face_eps   = 0.02;                 % [m] clearance
koz.bigM       = 1000;                  % [m]
koz.window     = [0.0, 1.0];            % active across whole horizon
koz.near_Tgate = inf;                   % |T|<=gate to activate

% Weights (match SCP spirit; can be tuned)
weights = struct();
weights.R   = 1e4;    % control effort (L2)
weights.w_tr= 1e-5;   % trust region (state) (L2)
weights.w_v = 1e2;    % virtual control (L2)

% Reference for trust region
[Xbar,Ubar] = local_qonly_seed(stm, N, RelInitState, RelFinalState, mission, dt);
if ~all(isfinite(Xbar(:))) || size(Xbar,2)~=N
    Xbar = zeros(6,N); for k=1:N, a=(k-1)/(N-1); Xbar(:,k)=(1-a)*RelInitState + a*RelFinalState; end
    Ubar = zeros(3,N-1);
end

% Build KOZ boxes
boxes = build_koz_aabbs(koz);
fprintf('KOZ AABBs generated: %d\n', size(boxes,3));

%% Build MIQP (Gurobi)
switch lower(solver)
    case 'gurobi'
        if exist('gurobi','file')~=2
            error('Gurobi MATLAB API not found on path. Set solver=''quadprog'' or add Gurobi to path.');
        end
        [model, idx] = build_miqp_gurobi(stm, N, RelInitState, RelFinalState, mission, weights, Xbar, boxes, koz);
        params_grb = struct(); params_grb.OutputFlag = 1; params_grb.MIPGap = 1e-3; params_grb.TimeLimit = 300; params_grb.Method = 1;
        result = gurobi(model, params_grb);
        if ~isfield(result,'x'); error('MIQP failed: status=%s', string(result.status)); end
        z = result.x; [X,U,V] = unpack(z, idx, N); solver_note = sprintf('Gurobi MIQP (gap=%s)', string(conditional_field(result,'mipgap','N/A')));
    otherwise % 'quadprog' relaxation
        [H,f,A,b,Aeq,beq,lb,ub,idx] = build_miqp_qprelax(stm, N, RelInitState, RelFinalState, mission, weights, Xbar, boxes, koz);
        opts = optimoptions('quadprog','Display','iter','MaxIterations',2000,'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8);
        [z, fval_qp, exitflag_qp] = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],opts); %#ok<ASGLU>
        [X,U,V] = unpack(z, idx, N); solver_note = 'quadprog QP-relaxation (binaries relaxed)';
end

%% Diagnostics and plots (match SCP_QP style)
dv = sum(vecnorm(U,2,1))*dt;
fprintf('Solver: %s, Cum dV: %.3f m/s\n', solver_note, dv);
local_visualize_like_scp(X, Xbar, U, time_vec, boxes);
try, print(gcf, fullfile(fileparts(mfilename('fullpath')), 'scp_koz_results_miqp.png'), '-dpng','-r150'); catch, end
try
    save(fullfile(fileparts(mfilename('fullpath')),'scp_miqp_last.mat'),'X','U','boxes','time_vec');
catch, end

%% ===== Helpers =====
function s = conditional_field(st, f, def)
    if isstruct(st) && isfield(st, f), s = string(st.(f)); else, s = string(def); end
end

function stm = local_computeSTM(params, chief, time_vec, RelInitState)
    B = zeros(6,3); B(2,1)=1; B(4,2)=1; B(6,3)=1; B_full=zeros(6,6); B_full(:,1:3)=B;
    initStruct.params = {params.Req, params.mu, params.J2, params.tol, params.safetyAltitude};
    initStruct.maneuverParams = {10, B_full};
    initStruct.timeParams.t0 = time_vec(1); initStruct.timeParams.dt = mean(diff(time_vec)); initStruct.timeParams.tf = time_vec(end);
    initStruct.initChiefDescription='Classical'; initStruct.initDeputyDescription='Cartesian'; initStruct.Elements = chief.elements; initStruct.RelInitState = RelInitState; initStruct.RelFinalState = [];
    stm = GimAlfriendSTM(initStruct); stm.makeDiscreteMatrices();
end

function [X_opt,U_opt] = local_qonly_seed(stm, N, x0, xf, mission, dt)
    nx=6; nu=3; nX=N*nx; nU=(N-1)*nu; nvar=nX+nU; Qx=diag([1,1e-3,1,1e-3,1,1e-3]);
    H=sparse(nvar,nvar); f=zeros(nvar,1); for k=1:N, ix=(k-1)*nx+(1:nx); H(ix,ix)=H(ix,ix)+2*Qx; end
    Aeq=sparse(nx+nx*(N-1)+nx, nvar); beq=zeros(size(Aeq,1),1); row=0;
    Aeq(1:nx,1:nx)=eye(nx); beq(1:nx)=x0; row=row+nx;
    for k=1:N-1, rows=row+(1:nx); Aeq(rows,k*nx+(1:nx))=eye(nx); Aeq(rows,(k-1)*nx+(1:nx))=-stm.Ak(:,:,k); Aeq(rows,nX+(k-1)*nu+(1:nu))=-stm.Bk(:,1:3,k); row=row+nx; end
    Aeq(row+(1:nx),(N-1)*nx+(1:nx))=eye(nx); beq(row+(1:nx))=xf;
    nI=2*(N-1)*nu; Ai=sparse(nI,nvar); bi=zeros(nI,1); r=0; for k=1:N-1, for j=1:nu, e=zeros(1,nu); e(j)=1; r=r+1; Ai(r,nX+(k-1)*nu+(1:nu))= e; bi(r)=mission.u_max; r=r+1; Ai(r,nX+(k-1)*nu+(1:nu))=-e; bi(r)=mission.u_max; end, end
    opts=optimoptions('quadprog','Display','off','MaxIterations',2000);
    z=quadprog(H,f,Ai,bi,Aeq,beq,[],[],[],opts);
    X_opt=reshape(z(1:nX),nx,N); U_opt=reshape(z(nX+1:end),nu,N-1);
end

function boxes = build_koz_aabbs(ko)
    boxes = zeros(2,3,0);
    if ~exist(ko.expr_dir,'dir')
        warning('Expr dir not found: %s', ko.expr_dir); return;
    end
    d = dir(fullfile(ko.expr_dir,'*_*_expr.txt'));
    take = min(numel(d), ko.max_expr);
    if take==0, return; end
    xv = linspace(ko.domain_lo(1), ko.domain_hi(1), ko.grid_res);
    yv = linspace(ko.domain_lo(2), ko.domain_hi(2), ko.grid_res);
    zv = linspace(ko.domain_lo(3), ko.domain_hi(3), ko.grid_res);
    [Xg,Yg,Zg] = meshgrid(xv,yv,zv);
    for i=1:take
        expr_path = fullfile(ko.expr_dir, d(i).name);
        Fh = make_F_from_expr(expr_path);
        if isempty(Fh), continue; end
        try
            V = arrayfun(@(x,y,z) Fh(x,y,z), Xg, Yg, Zg) - 1;
            S = isosurface(Xg, Yg, Zg, V, 0);
            if isempty(S.vertices), continue; end
            lo = min(S.vertices, [], 1); hi = max(S.vertices, [], 1);
            B = [lo(1) hi(1); lo(2) hi(2); lo(3) hi(3)]; % [x;y;z] rows
            boxes(:,:,end+1) = [B(1,:); B(2,:)]; %#ok<AGROW> % 2x3: low; high
        catch
            V = arrayfun(@(x,y,z) Fh(x,y,z), Xg, Yg, Zg) - 1;
            mask = V < 0; if ~any(mask(:)), continue; end
            [I,J,K] = ind2sub(size(V), find(mask));
            lo = [xv(min(I)) yv(min(J)) zv(min(K))];
            hi = [xv(max(I)) yv(max(J)) zv(max(K))];
            boxes(:,:,end+1) = [lo; hi]; %#ok<AGROW>
        end
    end
end

function Fh = make_F_from_expr(expr_path)
    try
        txt = strtrim(fileread(expr_path));
        txt = regexprep(txt, '\*\*', '.^'); % Python power to MATLAB
        fstr = sprintf('@(x,y,z) (%s)', txt);
        Fh = str2func(fstr);
    catch
        Fh = [];
    end
end

function [model, idx] = build_miqp_gurobi(stm, N, x0, xf, mission, weights, Xbar, boxes, ko)
    nx=6; nu=3; nv=6; Epos = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
    nX = N*nx; nU=(N-1)*nu; nV=(N-1)*nv;
    nBoxes = size(boxes,3);
    faces_per_box = 6;
    k0 = max(1, floor(ko.window(1)*(N-1))+1); k1 = min(N, ceil(ko.window(2)*(N-1))+1);
    nKwin = max(0, k1-k0+1);
    nZ = nKwin * nBoxes * faces_per_box;

    % variable order: [X, U, V, Zbin]
    nvars = nX + nU + nV + nZ;
    lb = -inf(nvars,1); ub = inf(nvars,1);
    vtype = repmat('C', nvars, 1);
    if nZ>0, vtype(end-nZ+1:end) = 'B'; end

    % indices
    idx.X = @(k) ((k-1)*nx+1):(k*nx);
    Xofs=0; Uofs = nX; Vofs = nX+nU; Zofs = nX+nU+nV;

    % bounds on U
    for k=1:N-1
        iu = Uofs+(k-1)*nu+(1:nu);
        lb(iu) = -mission.u_max; ub(iu) = mission.u_max;
    end

    % Objective Q (quadratic)
    Q = sparse(nvars, nvars);
    % control effort
    for k=1:N-1
        iu = Uofs+(k-1)*nu+(1:nu);
        Q(iu,iu) = Q(iu,iu) + 2*weights.R*eye(nu);
    end
    % virtual control
    for k=1:N-1
        iv = Vofs+(k-1)*nv+(1:nv);
        Q(iv,iv) = Q(iv,iv) + 2*weights.w_v*eye(nv);
    end
    % trust region on X around Xbar
    for k=1:N
        ix = Xofs+(k-1)*nx+(1:nx);
        Q(ix,ix) = Q(ix,ix) + 2*weights.w_tr*eye(nx);
    end
    obj = zeros(nvars,1);
    % linear term for trust region: -2*w_tr*Xbar' x  -> put in obj
    for k=1:N
        ix = Xofs+(k-1)*nx+(1:nx);
        obj(ix) = obj(ix) - 2*weights.w_tr*Xbar(:,k);
    end

    % Equalities: initial + dynamics + final
    nEq = nx + nx*(N-1) + nx;
    Aeq = sparse(nEq, nvars); rhs = zeros(nEq,1); row=0;
    % initial
    Aeq(row+(1:nx), idx.X(1)) = eye(nx); rhs(row+(1:nx)) = x0; row=row+nx;
    % dynamics
    for k=1:N-1
        rows = row+(1:nx);
        Aeq(rows, idx.X(k+1)) = eye(nx);
        Aeq(rows, idx.X(k))   = -stm.Ak(:,:,k);
        iu = Uofs+(k-1)*nu+(1:nu);
        Aeq(rows, iu)         = -stm.Bk(:,1:3,k);
        iv = Vofs+(k-1)*nv+(1:nv);
        Aeq(rows, iv)         = -eye(nx);
        row=row+nx;
    end
    % final
    Aeq(row+(1:nx), idx.X(N)) = eye(nx); rhs(row+(1:nx)) = xf;

    % Inequalities: AABB disjunction with big-M + sum(z)>=1
    Al = sparse(0, nvars); bl = zeros(0,1); sense = repmat(' ', 0, 1);
    zptr = 0;
    for k=k0:k1
        for i=1:nBoxes
            B = boxes(:,:,i); % 2x3: [low; high]
            lx=B(1,1); ux=B(2,1); ly=B(1,2); uy=B(2,2); lz=B(1,3); uz=B(2,3);
            pbar = Epos*Xbar(:,k); gate_on = (abs(pbar(2)) <= ko.near_Tgate);
            % z indices for this (k,i)
            z1 = Zofs+zptr+1; z2=z1+1; z3=z2+1; z4=z3+1; z5=z4+1; z6=z5+1;
            zidx = [z1 z2 z3 z4 z5 z6];
            if gate_on
                ix = idx.X(k);
                % x <= lx-eps + M*(1-z1)
                a = sparse(1, nvars); a(ix(1)) = 1;  a(z1) = ko.bigM; Al(end+1,:) = a; bl(end+1,1) = ko.bigM + (lx - ko.face_eps); sense(end+1,1) = '<';
                % x >= ux+eps - M*(1-z2)  -> -x <= -ux-eps + M*(1-z2)
                a = sparse(1, nvars); a(ix(1)) = -1; a(z2) = ko.bigM; Al(end+1,:) = a; bl(end+1,1) = ko.bigM - (ux + ko.face_eps); sense(end+1,1) = '<';
                % y <= ly-eps + M*(1-z3)
                a = sparse(1, nvars); a(ix(3)) = 1;  a(z3) = ko.bigM; Al(end+1,:) = a; bl(end+1,1) = ko.bigM + (ly - ko.face_eps); sense(end+1,1) = '<';
                % y >= uy+eps - M*(1-z4)
                a = sparse(1, nvars); a(ix(3)) = -1; a(z4) = ko.bigM; Al(end+1,:) = a; bl(end+1,1) = ko.bigM - (uy + ko.face_eps); sense(end+1,1) = '<';
                % z <= lz-eps + M*(1-z5)
                a = sparse(1, nvars); a(ix(5)) = 1;  a(z5) = ko.bigM; Al(end+1,:) = a; bl(end+1,1) = ko.bigM + (lz - ko.face_eps); sense(end+1,1) = '<';
                % z >= uz+eps - M*(1-z6)
                a = sparse(1, nvars); a(ix(5)) = -1; a(z6) = ko.bigM; Al(end+1,:) = a; bl(end+1,1) = ko.bigM - (uz + ko.face_eps); sense(end+1,1) = '<';
                % sum(z) >= 1  ->  -z1 - ... - z6 <= -1
                a = sparse(1, nvars); a(zidx) = -1; Al(end+1,:) = a; bl(end+1,1) = -1; sense(end+1,1) = '<';
            else
                % inactive gate: add no-ops with large RHS
                for t=1:faces_per_box+1, Al(end+1,1) = 0; bl(end+1,1) = 1e9; sense(end+1,1) = '<'; end %#ok<AGROW>
            end
            zptr = zptr + faces_per_box;
        end
    end

    % Assemble model
    model = struct();
    model.Q = Q; model.obj = obj; model.A = [Aeq; Al];
    model.rhs = [rhs; bl];
    model.sense = [repmat('=', size(Aeq,1),1); sense];
    model.lb = lb; model.ub = ub; model.vtype = vtype; model.modelsense = 'min';
end

function [H,f,A,b,Aeq,beq,lb,ub,idx] = build_miqp_qprelax(stm, N, x0, xf, mission, weights, Xbar, boxes, ko)
    % Same structure as build_miqp_gurobi, but returns quadprog inputs and relaxes binaries to [0,1]
    nx=6; nu=3; nv=6; Epos = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
    nX = N*nx; nU=(N-1)*nu; nV=(N-1)*nv;
    nBoxes = size(boxes,3); faces_per_box = 6;
    k0 = max(1, floor(ko.window(1)*(N-1))+1); k1 = min(N, ceil(ko.window(2)*(N-1))+1);
    nKwin = max(0, k1-k0+1); nZ = nKwin * nBoxes * faces_per_box;
    nvars = nX + nU + nV + nZ;
    lb = -inf(nvars,1); ub = inf(nvars,1);
    % bounds on U
    Uofs = nX; Vofs = nX+nU; Zofs = nX+nU+nV;
    for k=1:N-1, iu=Uofs+(k-1)*nu+(1:nu); lb(iu)=-mission.u_max; ub(iu)=mission.u_max; end
    % binaries relaxed to [0,1]
    if nZ>0, lb(Zofs+1:Zofs+nZ)=0; ub(Zofs+1:Zofs+nZ)=1; end
    % Quadratic H and linear f
    H = sparse(nvars,nvars); f = zeros(nvars,1);
    for k=1:N-1, iu=Uofs+(k-1)*nu+(1:nu); H(iu,iu)=H(iu,iu)+2*weights.R*eye(nu); end
    for k=1:N-1, iv=Vofs+(k-1)*nv+(1:nv); H(iv,iv)=H(iv,iv)+2*weights.w_v*eye(nv); end
    for k=1:N, ix=(k-1)*nx+(1:nx); H(ix,ix)=H(ix,ix)+2*weights.w_tr*eye(nx); f(ix)=f(ix)-2*weights.w_tr*Xbar(:,k); end
    % Equalities
    Aeq = sparse(nx + nx*(N-1) + nx, nvars); beq = zeros(size(Aeq,1),1); row=0;
    Aeq(row+(1:nx), 1:nx) = eye(nx); beq(row+(1:nx))=x0; row=row+nx;
    for k=1:N-1
        rows=row+(1:nx); Aeq(rows, (k)*nx+(1:nx))=eye(nx); Aeq(rows, (k-1)*nx+(1:nx))=-stm.Ak(:,:,k);
        iu = Uofs+(k-1)*nu+(1:nu); Aeq(rows, iu) = -stm.Bk(:,1:3,k);
        iv = Vofs+(k-1)*nv+(1:nv); Aeq(rows, iv) = -eye(nx);
        row=row+nx;
    end
    Aeq(row+(1:nx), (N-1)*nx+(1:nx))=eye(nx); beq(row+(1:nx))=xf;
    % Inequalities (<=)
    A = sparse(0,nvars); b = zeros(0,1);
    zptr = 0;
    for k=k0:k1
        for i=1:nBoxes
            B=boxes(:,:,i); lx=B(1,1); ux=B(2,1); ly=B(1,2); uy=B(2,2); lz=B(1,3); uz=B(2,3);
            pbar = Epos*Xbar(:,k); gate_on = (abs(pbar(2)) <= ko.near_Tgate);
            z1 = Zofs+zptr+1; z2=z1+1; z3=z2+1; z4=z3+1; z5=z4+1; z6=z5+1; zidx=[z1 z2 z3 z4 z5 z6];
            if gate_on
                ix=(k-1)*nx+(1:nx);
                a = sparse(1,nvars); a(ix(1))= 1; a(z1)=ko.bigM; A(end+1,:)=a; b(end+1,1)=ko.bigM + (lx - ko.face_eps); %#ok<AGROW>
                a = sparse(1,nvars); a(ix(1))=-1; a(z2)=ko.bigM; A(end+1,:)=a; b(end+1,1)=ko.bigM - (ux + ko.face_eps); %#ok<AGROW>
                a = sparse(1,nvars); a(ix(3))= 1; a(z3)=ko.bigM; A(end+1,:)=a; b(end+1,1)=ko.bigM + (ly - ko.face_eps); %#ok<AGROW>
                a = sparse(1,nvars); a(ix(3))=-1; a(z4)=ko.bigM; A(end+1,:)=a; b(end+1,1)=ko.bigM - (uy + ko.face_eps); %#ok<AGROW>
                a = sparse(1,nvars); a(ix(5))= 1; a(z5)=ko.bigM; A(end+1,:)=a; b(end+1,1)=ko.bigM + (lz - ko.face_eps); %#ok<AGROW>
                a = sparse(1,nvars); a(ix(5))=-1; a(z6)=ko.bigM; A(end+1,:)=a; b(end+1,1)=ko.bigM - (uz + ko.face_eps); %#ok<AGROW>
                a = sparse(1,nvars); a(zidx) = -1; A(end+1,:)=a; b(end+1,1) = -1; %#ok<AGROW>
            end
            zptr = zptr + faces_per_box;
        end
    end
    idx = struct(); idx.X=@(k)((k-1)*nx+1):(k*nx);
    lb = lb; ub = ub; % return
end

function [X,U,V] = unpack(z, idx, N)
    nx=6; nu=3; nv=6; nX=N*nx; nU=(N-1)*nu; nV=(N-1)*nv;
    X = reshape(z(1:nX), nx, N);
    U = reshape(z(nX+(1:nU)), nu, N-1);
    V = reshape(z(nX+nU+(1:nV)), nv, N-1);
end

function local_visualize_like_scp(X, Xbar, U, time_vec, boxes)
    figure('Position',[60 60 1200 700]);
    % 3D trajectory + KOZ
    subplot(2,2,1);
    hold on; grid on; box on;
    if nargin>=5 && ~isempty(boxes)
        for i=1:size(boxes,3)
            draw_box(boxes(:,:,i));
        end
    end
    if ~isempty(Xbar)
        h_init = plot3(Xbar(1,:)/1e3, Xbar(3,:)/1e3, Xbar(5,:)/1e3,'--','Color',[0.5 0.5 0.5],'LineWidth',1.5);
    else
        h_init = gobjects(0);
    end
    h_final = plot3(X(1,:)/1e3, X(3,:)/1e3, X(5,:)/1e3,'b-','LineWidth',2);
    h_start = plot3(X(1,1)/1e3, X(3,1)/1e3, X(5,1)/1e3,'go','MarkerFaceColor','g');
    h_goal  = plot3(X(1,end)/1e3, X(3,end)/1e3, X(5,end)/1e3,'ro','MarkerFaceColor','r');
    axis equal; view(45,30);
    xlim([-0.1 0.1]); ylim([-0.1 0.1]); zlim([-0.1 0.1]);
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]'); title('3D Trajectory (pre/post) + KOZ');
    if ~isempty(h_init)
        legend([h_init h_final h_start h_goal], {'Initial traj','Final traj','Start','Goal'}, 'Location','northwest');
    else
        legend([h_final h_start h_goal], {'Final traj','Start','Goal'}, 'Location','northwest');
    end
    % Control magnitude
    subplot(2,2,2);
    stairs(time_vec(1:end-1)/3600, vecnorm(U,2,1)*1e3,'k-','LineWidth',1.5); grid on;
    xlabel('Time [h]'); ylabel('|u| [mm/s^2]'); title('Control magnitude');
    % Range over time
    subplot(2,2,3);
    rel = sqrt(sum(X([1,3,5],:).^2,1));
    plot(time_vec/3600, rel,'m-','LineWidth',2); grid on; xlabel('Time [h]'); ylabel('Range [m]');
    % Cumulative dV
    subplot(2,2,4);
    dt = mean(diff(time_vec));
    dv_cum = [0 cumsum(vecnorm(U,2,1)*dt)];
    plot(time_vec/3600, dv_cum,'g-','LineWidth',2); grid on; xlabel('Time [h]'); ylabel('Cum. dV [m/s]');
end

function draw_box(B)
    % B: 2x3: [low; high] for x,y,z
    lx=B(1,1); ux=B(2,1); ly=B(1,2); uy=B(2,2); lz=B(1,3); uz=B(2,3);
    % 8 corners
    V = [lx ly lz;
         ux ly lz;
         ux uy lz;
         lx uy lz;
         lx ly uz;
         ux ly uz;
         ux uy uz;
         lx uy uz];
    V = V/1e3; % km for plotting
    % Faces (quads)
    F = [1 2 3 4;  % bottom z=lz
         5 6 7 8;  % top z=uz
         1 2 6 5;  % side x
         2 3 7 6;  % side y
         3 4 8 7;  % side x
         4 1 5 8]; % side y
    try
        patch('Faces',F,'Vertices',V,'FaceAlpha',0.05,'EdgeColor',[1 0 0], 'FaceColor',[1 0 0]);
    catch
        % Fallback: wireframe
        E = [1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
        hold on;
        for i=1:size(E,1)
            p1 = V(E(i,1),:); p2 = V(E(i,2),:);
            plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'Color',[1 0 0]);
        end
    end
end
