%% MILP-based KOZ avoidance with L1 objective
% - Dynamics: Gim–Alfriend STM (J2) in LVLH (RTN) coordinates
% - Objective: L1 control + L1 virtual control + optional trust region
% - Constraints:
%     x_{k+1} = A_k x_k + B_k u_k + v_k
%     |u_i| <= u_max (component-wise)
%     AABB KOZ (per implicit F file): outside-of-box via big-M disjunction
%       For each box B=[lx ux; ly uy; lz uz], enforce at least one face:
%         x<=lx-eps OR x>=ux+eps OR y<=ly-eps OR y>=uy+eps OR z<=lz-eps OR z>=uz+eps
%       Encoded with binaries and big-M.

clear; clc; close all;

%% Inputs
start_goal_csv = 'custom_states.csv'; % two rows: [x vx y vy z vz] (meters)
dt = 60;    % [s]
N  = 101;   % nodes
time_vec = linspace(0, dt*(N-1), N);

% Physical parameters
params = struct();
params.Req = 6378.137e3; params.mu = 3.986004415e14; params.J2 = 1082.63e-6;
params.tol = 1e-12; params.safetyAltitude = 100e3;

% Chief orbit (for STM only)
chief = struct();
chief.altitude = 700e3; chief.a = params.Req + chief.altitude;
chief.e = 0.001; chief.i = deg2rad(98.2); chief.RAAN = deg2rad(30);
chief.w = deg2rad(0); chief.M = deg2rad(60);
chief.elements = [chief.a, chief.e, chief.i, chief.RAAN, chief.w, chief.M];

% Mission
mission = struct('u_max', 0.1); % [m/s^2] component-wise bound

% MILP options
koz_opts = struct();
koz_opts.expr_dir   = 'C:\\Users\\98kim\\Desktop\\Acta-Astronautica\\Funcs_ISS_expr';
koz_opts.domain_lo  = [-100; -100; -100];   % [m]
koz_opts.domain_hi  = [ 100;  100;  100];   % [m]
koz_opts.grid_res   = 40;                   % grid per dimension for isosurface bbox
koz_opts.max_expr   = 12;                   % cap number of expr files
koz_opts.face_eps   = 0.02;                 % [m] clearance outside the box
koz_opts.bigM       = 1000;                 % big-M (meter scale ~100 -> M=1000 ok)
koz_opts.window     = [0.0, 1.0];           % active horizon fraction
koz_opts.near_Tgate = inf;                  % optional |T|<=gate activation

weights = struct('R',1e3, 'w_v',1e2, 'w_tr',1e-3); % L1 weights
use_trust_region = true;                     % include L1 trust region |x-xbar|

%% Load start/goal and STM
[x0, xf] = load_custom_states_milp(start_goal_csv);
fprintf('Start (m): %s\n', sprintf('% .3f ', x0));
fprintf('Goal  (m): %s\n', sprintf('% .3f ', xf));

stm = local_computeSTM(params, chief, time_vec, x0);

%% Initial guess (Q-only seed)
[Xbar,Ubar] = local_qonly_seed(stm, N, x0, xf, mission, dt);
if ~all(isfinite(Xbar(:))) || size(Xbar,2)~=N
    % Fall back to straight interpolation and zero U
    Xbar = zeros(6,N); for k=1:N, a=(k-1)/(N-1); Xbar(:,k)=(1-a)*x0 + a*xf; end; Ubar=zeros(3,N-1);
end

%% Build KOZ AABB list from implicit expr files
boxes = build_koz_aabbs(koz_opts);
fprintf('KOZ AABBs generated: %d\n', size(boxes,3));

%% Build MILP (intlinprog)
[c, A, b, Aeq, beq, lb, ub, intcon, idx] = build_milp(stm, N, x0, xf, mission, weights, Xbar, boxes, koz_opts, use_trust_region);

opts = optimoptions('intlinprog', 'Display','iter', 'MaxTime', 300, 'Heuristics','advanced', ...
    'LPMaxIter', 1e6, 'CutGeneration','intermediate');

fprintf('\nSolving MILP ...\n');
[z, fval, exitflag, output] = intlinprog(c, intcon, A, b, Aeq, beq, lb, ub, opts);

if exitflag <= 0
    fprintf('MILP failed: flag=%d\n', exitflag);
end

[X,U,V,Tu,Tv,Dx] = unpack_milp(z, idx, N); %#ok<ASGLU>

%% Diagnostics and plots
dv = sum(vecnorm(U,2,1))*dt;
fprintf('Objective(L1): %.6e, Cum dV: %.3f m/s, exit=%d, iters=%d\n', fval, dv, exitflag, getfield(output,'relativegap',NaN)); %#ok<GFLD>

figure('Position',[60 60 1200 700]);
subplot(2,2,1); hold on; grid on; box on;
% Draw boxes
for i=1:size(boxes,3)
    draw_box(boxes(:,:,i));
end
plot3(Xbar(1,:)/1e3, Xbar(3,:)/1e3, Xbar(5,:)/1e3,'--','Color',[0.6 0.6 0.6],'LineWidth',1.0);
plot3(X(1,:)/1e3, X(3,:)/1e3, X(5,:)/1e3,'b-','LineWidth',2);
plot3(X(1,1)/1e3, X(3,1)/1e3, X(5,1)/1e3,'go','MarkerFaceColor','g');
plot3(X(1,end)/1e3, X(3,end)/1e3, X(5,end)/1e3,'ro','MarkerFaceColor','r');
axis equal; view(45,30);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]'); title('MILP Trajectory with AABB KOZ');

subplot(2,2,2);
stairs(time_vec(1:end-1)/3600, vecnorm(U,2,1)*1e3,'k-','LineWidth',1.5); grid on;
xlabel('Time [h]'); ylabel('|u| [mm/s^2]'); title('Control magnitude');

subplot(2,2,3);
rel = sqrt(sum(X([1,3,5],:).^2,1));
plot(time_vec/3600, rel,'m-','LineWidth',2); grid on; xlabel('Time [h]'); ylabel('Range [m]');

subplot(2,2,4);
dv_cum = [0 cumsum(vecnorm(U,2,1)*dt)];
plot(time_vec/3600, dv_cum,'g-','LineWidth',2); grid on; xlabel('Time [h]'); ylabel('Cum. dV [m/s]');

try
    save(fullfile(fileparts(mfilename('fullpath')),'scp_milp_last.mat'),'X','U','boxes','time_vec');
catch, end

%% ===== Helpers =====
function [x0, xf] = load_custom_states_milp(csvp)
    if exist(csvp,'file')
        Xio = readmatrix(csvp); x0 = Xio(1,:).'; xf = Xio(2,:).';
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
            B = [lo(1) hi(1); lo(2) hi(2); lo(3) hi(3)];
            boxes(:,:,end+1) = B.'; %#ok<AGROW>
        catch
            % fallback: coarse bbox via threshold
            V = arrayfun(@(x,y,z) Fh(x,y,z), Xg, Yg, Zg) - 1;
            mask = V < 0; if ~any(mask(:)), continue; end
            [I,J,K] = ind2sub(size(V), find(mask));
            lo = [xv(min(I)) yv(min(J)) zv(min(K))];
            hi = [xv(max(I)) yv(max(J)) zv(max(K))];
            B = [lo(1) hi(1); lo(2) hi(2); lo(3) hi(3)];
            boxes(:,:,end+1) = B.'; %#ok<AGROW>
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

function [c, A, b, Aeq, beq, lb, ub, intcon, idx] = build_milp(stm, N, x0, xf, mission, weights, Xbar, boxes, ko, use_tr)
    nx=6; nu=3; nv=6; nX=N*nx; nU=(N-1)*nu; nV=(N-1)*nv;
    nTU=(N-1)*nu; nTV=(N-1)*nv; nDX = N*nx;
    nBoxes = size(boxes,3);
    faces_per_box = 6; % x<=lx, x>=ux, y<=ly, y>=uy, z<=lz, z>=uz
    % active window
    k0 = max(1, floor(ko.window(1)*(N-1))+1); k1 = min(N, ceil(ko.window(2)*(N-1))+1);
    % count binaries (only in window)
    nZ = (k1-k0+1) * nBoxes * faces_per_box;

    % variable ordering
    % z = [X | U | V | TU | TV | DX | Zbin]
    nvar = nX+nU+nV+nTU+nTV+nDX + nZ;
    c    = zeros(nvar,1);
    lb   = -inf(nvar,1); ub =  inf(nvar,1);

    % indices
    idx.X  = @(k) ((k-1)*nx+1):(k*nx);
    offset = 0; Xofs=offset; offset=offset+nX;
    idxU = @(k) (Xofs+nX + (k-1)*nu+1):(Xofs+nX + k*nu); % helper not used
    idx.U  = @(k) (nX+(k-1)*nu+1):(nX+k*nu);
    Uofs = nX; Vofs = nX+nU; TUofs = nX+nU+nV; TVofs = TUofs+nTU; DXofs = TVofs+nTV; Zofs = DXofs+nDX;
    % L1 objective weights
    for k=1:N-1
        iu = (TUofs+(k-1)*nu+1):(TUofs+k*nu);
        c(iu) = c(iu) + weights.R*ones(nu,1);
        iv = (TVofs+(k-1)*nv+1):(TVofs+k*nv);
        c(iv) = c(iv) + weights.w_v*ones(nv,1);
    end
    if use_tr
        for k=1:N
            ix = (DXofs+(k-1)*nx+1):(DXofs+k*nx);
            c(ix) = c(ix) + weights.w_tr*ones(nx,1);
        end
    end

    % bounds
    lb(nX+1:nX+nU) = -mission.u_max; ub(nX+1:nX+nU) = mission.u_max; % u bounds
    lb(TUofs+1:TUofs+nTU) = 0; lb(TVofs+1:TVofs+nTV) = 0; lb(DXofs+1:DXofs+nDX) = 0;
    % binaries
    intcon = (Zofs+1):(Zofs+nZ);
    lb(intcon) = 0; ub(intcon) = 1;

    % Equality: initial + dynamics + final
    nEq = nx + nx*(N-1) + nx;
    Aeq = sparse(nEq, nvar); beq = zeros(nEq,1); row=0;
    % initial
    Aeq(row+(1:nx), idx.X(1)) = eye(nx); beq(row+(1:nx)) = x0; row=row+nx;
    % dynamics
    for k=1:N-1
        rows = row+(1:nx); Aeq(rows, idx.X(k+1))=eye(nx);
        Aeq(rows, idx.X(k)) = -stm.Ak(:,:,k);
        Aeq(rows, (nX+(k-1)*nu+1):(nX+k*nu)) = -stm.Bk(:,1:3,k);
        Aeq(rows, (Vofs+(k-1)*nv+1):(Vofs+k*nv)) = -eye(nx);
        row=row+nx;
    end
    % final
    Aeq(row+(1:nx), idx.X(N)) = eye(nx); beq(row+(1:nx)) = xf;

    % Inequalities
    % 1) TU >= +/- U  and TV >= +/- V,  DX >= |X - Xbar|
    % count rows
    nI = 2*nTU + 2*nTV + (use_tr*2*nDX);
    % 2) AABB KOZ disjunction rows: for each (k,i), 6 big-M + 1 sum(z)>=1
    nKwin = (k1-k0+1); nPairs = nKwin * nBoxes; nI = nI + nPairs*(faces_per_box + 1);
    Ai = sparse(nI, nvar); bi = zeros(nI,1); r=0;
    % TU/TV
    for k=1:N-1
        iu = (nX+(k-1)*nu+1):(nX+k*nu); tu = (TUofs+(k-1)*nu+1):(TUofs+k*nu);
        % TU >= U
        r=r+1; Ai(r,tu) = -eye(nu,nu); Ai(r,iu) = eye(nu,nu); bi(r)=0; r=r+nu-1;
        % TU >= -U
        r=r+1; Ai(r,tu) = -eye(nu,nu); Ai(r,iu) = -eye(nu,nu); bi(r)=0; r=r+nu-1;
        iv = (Vofs+(k-1)*nv+1):(Vofs+k*nv); tv = (TVofs+(k-1)*nv+1):(TVofs+k*nv);
        % TV >= V
        r=r+1; Ai(r,tv) = -eye(nv,nv); Ai(r,iv) = eye(nv,nv); bi(r)=0; r=r+nv-1;
        % TV >= -V
        r=r+1; Ai(r,tv) = -eye(nv,nv); Ai(r,iv) = -eye(nv,nv); bi(r)=0; r=r+nv-1;
    end
    % DX >= |X - Xbar|
    if use_tr
        for k=1:N
            ix = idx.X(k); dx = (DXofs+(k-1)*nx+1):(DXofs+k*nx);
            r=r+1; Ai(r,dx) = -eye(nx,nx); Ai(r,ix) = eye(nx,nx); bi(r) = Xbar(:,k); r=r+nx-1;
            r=r+1; Ai(r,dx) = -eye(nx,nx); Ai(r,ix) = -eye(nx,nx); bi(r) = -Xbar(:,k); r=r+nx-1;
        end
    end
    % 3) AABB KOZ with big-M
    % position extractor
    Epos = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
    zrow_start = r;
    zptr = 0; % for mapping (k,i,face)->global z index
    for k=k0:k1
        xk = idx.X(k);
        for i=1:nBoxes
            B = boxes(:,:,i); lx=B(1,1); ux=B(1,2); ly=B(2,1); uy=B(2,2); lz=B(3,1); uz=B(3,2);
            % gate by |T|<=near_Tgate
            pbar = Epos*Xbar(:,k);
            gate_on = (abs(pbar(2)) <= ko.near_Tgate);
            if gate_on
                % build 6 big-M inequalities
                % 1) x <= lx - eps + M*(1 - z1)
                z1 = Zofs + zptr + 1; z2 = z1+1; z3=z2+1; z4=z3+1; z5=z4+1; z6=z5+1;
                % x
                r=r+1; Ai(r, xk(1)) =  1; Ai(r, z1) =  ko.bigM; bi(r) =  ko.bigM + (lx - ko.face_eps);
                r=r+1; Ai(r, xk(1)) = -1; Ai(r, z2) =  ko.bigM; bi(r) =  ko.bigM - (ux + ko.face_eps);
                % y
                r=r+1; Ai(r, xk(3)) =  1; Ai(r, z3) =  ko.bigM; bi(r) =  ko.bigM + (ly - ko.face_eps);
                r=r+1; Ai(r, xk(3)) = -1; Ai(r, z4) =  ko.bigM; bi(r) =  ko.bigM - (uy + ko.face_eps);
                % z
                r=r+1; Ai(r, xk(5)) =  1; Ai(r, z5) =  ko.bigM; bi(r) =  ko.bigM + (lz - ko.face_eps);
                r=r+1; Ai(r, xk(5)) = -1; Ai(r, z6) =  ko.bigM; bi(r) =  ko.bigM - (uz + ko.face_eps);
                % sum(z) >= 1  -> -z1 - ... - z6 <= -1
                r=r+1; Ai(r, z1:z6) = -1; bi(r) = -1;
            else
                % inactive window: no-op with large RHS
                r=r+faces_per_box+1;
            end
            zptr = zptr + faces_per_box;
        end
    end
    % shrink A, b if we overshot
    if r < nI
        Ai = Ai(1:r, :); bi = bi(1:r);
    end
    A = Ai; b = bi;
end

function [X,U,V,Tu,Tv,Dx] = unpack_milp(z, idx, N)
    nx=6; nu=3; nv=6; nX=N*nx; nU=(N-1)*nu; nV=(N-1)*nv; nTU=(N-1)*nu; nTV=(N-1)*nv; nDX=N*nx;
    X = reshape(z(1:nX), nx, N);
    o = nX; U = reshape(z(o+(1:nU)), nu, N-1);
    o = o + nU; V = reshape(z(o+(1:nV)), nv, N-1);
    o = o + nV; Tu = reshape(z(o+(1:nTU)), nu, N-1);
    o = o + nTU; Tv = reshape(z(o+(1:nTV)), nv, N-1);
    o = o + nTV; Dx = reshape(z(o+(1:nDX)), nx, N);
end

function draw_box(B)
    % B: 2x3: [lx ux; ly uy; lz uz]'
    lx=B(1,1); ux=B(1,2); ly=B(2,1); uy=B(2,2); lz=B(3,1); uz=B(3,2);
    [Xc,Yc,Zc] = ndgrid([lx ux],[ly uy],[lz uz]);
    P = [Xc(:), Yc(:), Zc(:)];
    K = convhull(P);
    trisurf(K, P(:,1)/1e3, P(:,2)/1e3, P(:,3)/1e3, 'FaceAlpha',0.05, 'EdgeColor',[1 0 0], 'FaceColor',[1 0 0]);
end

