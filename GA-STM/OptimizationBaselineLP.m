%% 연료 최적 랑데부 - LP with Real Thruster Model using Gurobi
% 지구 텍스처, 개선된 케플러 방정식, Deputy 설정 모드, 실제 추력기 모델
% ECI 변환 오류 수정 버전
clear; clc; close all;

%% Gurobi 확인
if ~exist('gurobi', 'file')
    error('Gurobi가 설치되어 있지 않습니다! Gurobi를 설치하고 MATLAB 경로에 추가하세요.');
end
fprintf('Gurobi가 정상적으로 감지되었습니다.\n\n');

%% 1. 시스템 파라미터 설정
params = struct();
params.Req = 6378.137e3;           % 지구 반지름 [m]
params.mu = 3.986004415e14;        % 중력 상수 [m^3/s^2]
params.J2 = 1082.63e-6;            % J2 섭동 계수
params.tol = 1e-12;                % 케플러 방정식 허용 오차
params.safetyAltitude = 100e3;     % 안전 고도 [m]

% 주위성 궤도 요소
chief = struct();
chief.altitude = 550e3;            % 고도 [m]
chief.a = params.Req + chief.altitude;
chief.e = 0.001;                   % 이심률 (거의 원형 궤도)
chief.i = deg2rad(45.0);           % 궤도 경사각 [rad]
chief.RAAN = deg2rad(100.0);       % 승교점 적경 [rad]
chief.w = deg2rad(0);              % 근지점 편각 [rad] (원형 궤도이므로 0으로 설정)
chief.M = deg2rad(0);              % 평균 근점 이각 [rad] (임무 시작 시 기준점으로 설정)
chief.elements = [chief.a, chief.e, chief.i, chief.RAAN, chief.w, chief.M];
chief.n = sqrt(params.mu/chief.a^3);
chief.T = 2*pi/chief.n;

%% 2. 미션 파라미터
mission = struct();
mission.r0 = 50000;                % 초기 거리 [m]
mission.rf = 20;                   % 최종 거리 [m]

% 각 축별 독립적인 추력기 최대 추력 (물리적으로 정확한 모델)
mission.u_max_x = 1;            % X축 추력기 최대 추력 [m/s^2]
mission.u_max_y = 1;            % Y축 추력기 최대 추력 [m/s^2]
mission.u_max_z = 1;            % Z축 추력기 최대 추력 [m/s^2]

mission.constraint_type = 'hard';   % 'hard' or 'soft'
mission.final_state_tol = 1e-3;    % Hard constraint 허용 오차 [m]

% Deputy 초기화 모드 선택
deputy_init_mode = 'OrbitalElements';  % 'Cartesian' 또는 'OrbitalElements'

if strcmp(deputy_init_mode, 'Cartesian')
    % 기존 Cartesian 모드
    [RelInitState, RelFinalState] = setupInitialFinalStates_enhanced(mission.r0, mission.rf, 'Cartesian', [], chief, params);
    
elseif strcmp(deputy_init_mode, 'OrbitalElements')
    % 궤도 요소 기반 모드
    deputy_elements = struct();
    
    % Deputy의 초기 궤도 요소 차이 설정
    deputy_elements.initial.da = 50000;         % Semi-major axis 차이 [m] (+50km, 더 높은 에너지의 궤도)
    deputy_elements.initial.de = 0.019;         % Eccentricity 차이 (0.02 - 0.001, 목표가 더 뚜렷한 타원 궤도)
    deputy_elements.initial.di = deg2rad(0.2);      % Inclination 차이 [rad] (+0.2°, 궤도면 기울기 차이)
    deputy_elements.initial.dRAAN = deg2rad(5.0);     % RAAN 차이 [rad] (+5.0°, J2 섭동으로 인한 큰 궤도면 회전 차이)
    deputy_elements.initial.dw = deg2rad(1.0);      % Argument of perigee 차이 [rad] (궤도 타원의 방향 차이)
    deputy_elements.initial.dM = deg2rad(40.0);     % Mean anomaly 차이 [rad] (+40°, 임무 시작 시 위상차)
    
    [RelInitState, RelFinalState] = setupInitialFinalStates_enhanced(mission.r0, mission.rf, ...
                                                                    'OrbitalElements', deputy_elements, chief, params);
    
    fprintf('\nDeputy 초기 궤도 요소 차이:\n');
    fprintf('  Δa = %.1f m\n', deputy_elements.initial.da);
    fprintf('  Δe = %.6f\n', deputy_elements.initial.de);
    fprintf('  Δi = %.4f°\n', rad2deg(deputy_elements.initial.di));
    fprintf('  ΔRAAN = %.4f°\n', rad2deg(deputy_elements.initial.dRAAN));
    fprintf('  Δω = %.4f°\n', rad2deg(deputy_elements.initial.dw));
    fprintf('  ΔM = %.4f°\n', rad2deg(deputy_elements.initial.dM));
end

% 시간 설정
time_opts = struct();
time_opts.total_orbits = 100;
time_opts.dt = 120;
total_time = time_opts.total_orbits * chief.T;
time_vec = 0:time_opts.dt:total_time;
N = length(time_vec);

fprintf('\n=== 연료 최적 랑데부 (LP with Real Thruster Model) ===\n');
fprintf('Deputy 초기화 모드: %s\n', deputy_init_mode);
fprintf('초기 거리: %.1f km\n', mission.r0/1000);
fprintf('최종 거리: %.1f m (Hard Constraint)\n', mission.rf);
fprintf('최대 추력 (X/Y/Z): %.1f/%.1f/%.1f mm/s²\n', ...
        mission.u_max_x*1000, mission.u_max_y*1000, mission.u_max_z*1000);
fprintf('시간 스텝: %d\n', N);

%% 3. STM 계산
fprintf('\nSTM 계산 중...\n');
tic;
stm = computeSTM(params, chief, time_vec, RelInitState);
time_stm = toc;
fprintf('  완료: %.3f 초\n\n', time_stm);

%% 4. 최적화 문제 설정 - LP with Real Thruster Model
fprintf('최적화 문제 구성 중 (LP with Real Thruster Model)...\n');

% 최적화 방법 선택
optimization_method = 'LP';  % 'LP' 또는 'QP_comparison'

if strcmp(optimization_method, 'LP')
    % LP로 실제 추력기 모델 최적화
    lambda_values = [0];  % LP에서는 regularization 불필요
    results_all = cell(length(lambda_values), 1);
    
    for idx = 1:length(lambda_values)
        fprintf('\n--- 실제 추력기 모델 (LP) ---\n');
        
        % Gurobi로 LP 문제 해결
        tic;
        [X_opt, U_opt, exitflag, solve_info] = solveLP_RealThrusters(stm, N, RelInitState, ...
                                                             RelFinalState, mission, time_opts);
        time_solve = toc;
        
        % 결과 계산
        total_dV = computeTotalDeltaV_BoxConstraints(U_opt, time_opts.dt);
        final_error = norm(X_opt(:,end) - RelFinalState);
        final_pos_error = norm(X_opt([1,3,5],end) - RelFinalState([1,3,5]));
        final_vel_error = norm(X_opt([2,4,6],end) - RelFinalState([2,4,6]));
        
        fprintf('  최적화 완료: %.3f 초 (exitflag: %d)\n', time_solve, exitflag);
        fprintf('  Gurobi 실제 해결 시간: %.3f 초\n', solve_info.runtime);
        fprintf('  총 ΔV: %.3f m/s\n', total_dV);
        fprintf('  최종 위치 오차: %.2e m\n', final_pos_error);
        fprintf('  최종 속도 오차: %.2e m/s\n', final_vel_error);
        
        % 결과 저장
        results_all{idx} = struct('method', 'LP', 'X_opt', X_opt, 'U_opt', U_opt, ...
                                 'total_dV', total_dV, ...
                                 'final_error', final_error, 'exitflag', exitflag);
    end
    
elseif strcmp(optimization_method, 'QP_comparison')
    % QP와 LP 비교
    fprintf('\n=== LP vs QP 직접 비교 ===\n');
    
    % LP 해결 (실제 추력기 모델)
    fprintf('\n1. LP 해결 (실제 추력기 모델)...\n');
    [X_lp, U_lp, exitflag_lp, solve_info_lp] = solveLP_RealThrusters(stm, N, RelInitState, ...
                                                                  RelFinalState, mission, time_opts);
    total_dV_lp = computeTotalDeltaV_BoxConstraints(U_lp, time_opts.dt);
    
    % QP 해결
    fprintf('\n2. QP 해결...\n');
    weights.R = 1;
    [X_qp, U_qp, exitflag_qp] = solveQP_Gurobi(stm, N, RelInitState, RelFinalState, mission, weights);
    total_dV_qp = computeTotalDeltaV(U_qp, time_opts.dt);
    
    fprintf('\n총 ΔV - LP: %.3f m/s, QP: %.3f m/s (차이: %.1f%%)\n', ...
            total_dV_lp, total_dV_qp, (total_dV_qp - total_dV_lp) / total_dV_lp * 100);
    
    % 비교용 results_all 생성
    results_all = cell(2, 1);
    results_all{1} = struct('method', 'LP', 'X_opt', X_lp, 'U_opt', U_lp, ...
                           'total_dV', total_dV_lp, 'final_error', norm(X_lp(:,end) - RelFinalState));
    results_all{2} = struct('method', 'QP', 'X_opt', X_qp, 'U_opt', U_qp, ...
                           'total_dV', total_dV_qp, 'final_error', norm(X_qp(:,end) - RelFinalState));
end

%% 5. 결과 시각화 (LVLH Frame)
visualizeComparison_SOCP(results_all, time_vec, mission, chief);

%% 6. ECI 좌표계 시각화 (향상된 버전)
fprintf('\n=== ECI 좌표계 시각화 ===\n');
visualizeECITrajectories_enhanced(results_all, time_vec, mission, chief, params);

%% 보조 함수들

function [r_lvlh, v_lvlh] = eci_to_lvlh_relative(r_rel_eci, v_rel_eci, r_chief_eci, v_chief_eci)
    % ECI 기준 상대 상태를 LVLH 기준으로 변환하는 함수

    % 1. Chief의 위치/속도로부터 LVLH 프레임의 회전 행렬 계산
    r_norm = norm(r_chief_eci);
    R_hat = r_chief_eci / r_norm; % Radial (x-axis)
    
    h_vec = cross(r_chief_eci, v_chief_eci);
    h_norm = norm(h_vec);
    W_hat = h_vec / h_norm; % Cross-track (z-axis)
    
    S_hat = cross(W_hat, R_hat); % Along-track (y-axis)
    
    % ECI to LVLH 회전 행렬 (LVLH to ECI 행렬의 전치 행렬)
    C_eci2lvlh = [R_hat'; S_hat'; W_hat'];

    % 2. 위치 벡터 변환
    r_lvlh = C_eci2lvlh * r_rel_eci;

    % 3. 속도 벡터 변환 (코리올리 효과 포함)
    omega = h_vec / r_norm^2; % 궤도 각속도 벡터
    % 수정: 코리올리 효과의 올바른 계산
    v_lvlh = C_eci2lvlh * (v_rel_eci - cross(omega, r_rel_eci));
end

function [RelInitState, RelFinalState] = setupInitialFinalStates_enhanced(r0, rf, mode, deputy_elements, chief, params)
    % 향상된 초기/최종 상태 설정 (Cartesian 또는 정확한 Orbital Elements 모드)
    % chief와 params를 인자로 받도록 수정
    
    if strcmp(mode, 'Cartesian')
        % 기존 Cartesian 모드
        dx0 = r0 * 0.7071;
        dy0 = r0 * 0.7071;
        dz0 = 10000;
        RelInitState = [dx0; 0; dy0; 0; dz0; 0];
        
    elseif strcmp(mode, 'OrbitalElements')
        % === 정확한 궤도 요소 기반 모드 ===
        fprintf('\n정확한 궤도 요소 변환을 통해 초기 상태 계산 중...\n');
        
        % Deputy 궤도 요소 = Chief 요소 + 차이
        deputy.a = chief.a + deputy_elements.initial.da;
        deputy.e = chief.e + deputy_elements.initial.de;
        deputy.i = chief.i + deputy_elements.initial.di;
        deputy.RAAN = chief.RAAN + deputy_elements.initial.dRAAN;
        deputy.w = chief.w + deputy_elements.initial.dw;
        deputy.M = chief.M + deputy_elements.initial.dM;
        deputy.elements = [deputy.a, deputy.e, deputy.i, deputy.RAAN, deputy.w, deputy.M];

        % 2. t=0에서 각 위성의 ECI 상태 벡터 계산
        % 기존에 만들어진 keplerian2ECI_improved 함수 활용
        [r_chief_eci_init, v_chief_eci_init] = keplerian2ECI_improved(chief.elements, params.mu, 0);
        [r_deputy_eci_init, v_deputy_eci_init] = keplerian2ECI_improved(deputy.elements, params.mu, 0);

        % 3. ECI 기준 상대 상태 계산
        r_rel_eci = r_deputy_eci_init - r_chief_eci_init;
        v_rel_eci = v_deputy_eci_init - v_chief_eci_init;
        
        % 4. ECI 상대 상태를 LVLH 상대 상태로 변환
        [r_rel_lvlh, v_rel_lvlh] = eci_to_lvlh_relative(r_rel_eci, v_rel_eci, r_chief_eci_init, v_chief_eci_init);
        
        % 최종적으로 시뮬레이션에 사용될 초기 상태 벡터 (LVLH 기준)
        RelInitState = [r_rel_lvlh(1); v_rel_lvlh(1); r_rel_lvlh(2); v_rel_lvlh(2); r_rel_lvlh(3); v_rel_lvlh(3)];
        
        fprintf('  초기 상대 위치 (LVLH): [x: %.1f, y: %.1f, z: %.1f] m\n', RelInitState(1), RelInitState(3), RelInitState(5));
        fprintf('  초기 상대 속도 (LVLH): [vx: %.3f, vy: %.3f, vz: %.3f] m/s\n', RelInitState(2), RelInitState(4), RelInitState(6));
    else
        error('Unknown mode: %s. Use "Cartesian" or "OrbitalElements"', mode);
    end
    
    % 최종 목표 상태는 동일
    dxf = rf * 0.7071;
    dyf = rf * 0.7071;
    dzf = 0;
    RelFinalState = [dxf; 0; dyf; 0; dzf; 0];
end

function stm = computeSTM(params, chief, time_vec, RelInitState)
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

function [X_opt, U_opt, exitflag, solve_info] = solveLP_RealThrusters(stm, N, RelInitState, ...
                                                          RelFinalState, mission, time_opts)
    % 실제 추력기를 모델링한 LP 문제 (박스 제약)
    
    n_states = 6;
    n_controls = 3;
    n_x = N * n_states;
    n_u = (N-1) * n_controls;
    
    % L1 norm을 위한 보조 변수 추가
    n_aux = (N-1) * n_controls;  % |u_x|, |u_y|, |u_z| 변수들
    n_vars = n_x + n_u + n_aux;
    
    % Gurobi LP 모델
    model = struct();
    model.modelsense = 'min';
    
    %% 목적 함수: sum(|u_x| + |u_y| + |u_z|) * dt
    model.obj = zeros(n_vars, 1);
    model.obj(n_x + n_u + 1 : end) = time_opts.dt;  % aux 변수들의 합
    
    %% 변수 경계
    model.lb = -inf(n_vars, 1);
    model.ub = inf(n_vars, 1);
    
    % 각 축별 추력 제약 (박스 제약)
    for k = 1:N-1
        idx_ux = n_x + (k-1)*n_controls + 1;
        idx_uy = n_x + (k-1)*n_controls + 2;
        idx_uz = n_x + (k-1)*n_controls + 3;
        
        model.lb(idx_ux) = -mission.u_max_x;
        model.ub(idx_ux) = mission.u_max_x;
        
        model.lb(idx_uy) = -mission.u_max_y;
        model.ub(idx_uy) = mission.u_max_y;
        
        model.lb(idx_uz) = -mission.u_max_z;
        model.ub(idx_uz) = mission.u_max_z;
    end
    
    % aux 변수는 양수
    model.lb(n_x + n_u + 1 : end) = 0;
    
    %% 제약 조건
    % 1. 동역학 제약 (등식)
    [Aeq, beq] = setupEqualityConstraints_LP(stm, N, RelInitState, RelFinalState, n_x + n_u);
    
    % Aeq를 확장하여 aux 변수 포함
    Aeq_extended = [Aeq, sparse(size(Aeq,1), n_aux)];
    
    % 2. |u_i| <= aux_i 제약 (부등식)
    n_abs_constraints = 2 * (N-1) * n_controls;
    A_abs = sparse(n_abs_constraints, n_vars);
    b_abs = zeros(n_abs_constraints, 1);
    
    row = 0;
    for k = 1:N-1
        for j = 1:n_controls
            idx_u = n_x + (k-1)*n_controls + j;
            idx_aux = n_x + n_u + (k-1)*n_controls + j;
            
            % u_j <= aux_j
            row = row + 1;
            A_abs(row, idx_u) = 1;
            A_abs(row, idx_aux) = -1;
            b_abs(row) = 0;
            
            % -u_j <= aux_j
            row = row + 1;
            A_abs(row, idx_u) = -1;
            A_abs(row, idx_aux) = -1;
            b_abs(row) = 0;
        end
    end
    
    % 모든 제약 결합
    model.A = [Aeq_extended; A_abs];
    model.rhs = [beq; b_abs];
    model.sense = [repmat('=', length(beq), 1); 
                   repmat('<', n_abs_constraints, 1)];
    
    %% Gurobi 파라미터 (LP는 매우 빠름!)
    params = struct();
    params.OutputFlag = 1;
    params.Method = 1;              % Dual simplex (LP에 최적)
    params.Presolve = 2;            % Aggressive presolve
    params.FeasibilityTol = 1e-6;   % 적절한 허용 오차
    params.OptimalityTol = 1e-6;
    params.Threads = 0;             % 모든 코어 사용
    
    %% 문제 해결
    fprintf('  Gurobi로 LP 해결 중 (N=%d, 변수=%d, 제약=%d)...\n', ...
            N, n_vars, size(model.A, 1));
    
    result = gurobi(model, params);
    
    % 결과 정보 저장
    solve_info = struct();
    solve_info.status = result.status;
    solve_info.runtime = result.runtime;
    if isfield(result, 'itercount')
        solve_info.iterations = result.itercount;
    end
    if isfield(result, 'objval')
        solve_info.objval = result.objval;
    end
    
    %% 결과 처리
    if strcmp(result.status, 'OPTIMAL')
        exitflag = 1;
        sol = result.x;
        
        % 해 추출
        X_opt = reshape(sol(1:n_x), n_states, N);
        U_opt = reshape(sol(n_x+1:n_x+n_u), n_controls, N-1);
        
    elseif strcmp(result.status, 'SUBOPTIMAL')
        exitflag = 2;
        warning('Gurobi가 suboptimal 해를 반환했습니다.');
        sol = result.x;
        
        X_opt = reshape(sol(1:n_x), n_states, N);
        U_opt = reshape(sol(n_x+1:n_x+n_u), n_controls, N-1);
        
    else
        exitflag = 0;
        warning('Gurobi 해결 실패: %s', result.status);
        X_opt = zeros(n_states, N);
        U_opt = zeros(n_controls, N-1);
    end
end

function [X_opt, U_opt, exitflag] = solveQP_Gurobi(stm, N, RelInitState, RelFinalState, mission, weights)
    % Gurobi로 QP 문제 해결 (비교용)
    
    n_states = 6;
    n_controls = 3;
    n_vars = N*n_states + (N-1)*n_controls;
    
    % Gurobi 모델
    model = struct();
    model.modelsense = 'min';
    
    % 이차 목적 함수: u'*R*u
    model.Q = sparse(n_vars, n_vars);
    model.obj = zeros(n_vars, 1);
    
    R = weights.R * eye(3);
    for k = 1:N-1
        idx = N*n_states + (k-1)*n_controls + (1:n_controls);
        model.Q(idx, idx) = 2*R;
    end
    
    % 변수 경계
    model.lb = -inf(n_vars, 1);
    model.ub = inf(n_vars, 1);
    
    % 제어 입력 제약 (박스 제약)
    for k = 1:N-1
        idx_ux = N*n_states + (k-1)*n_controls + 1;
        idx_uy = N*n_states + (k-1)*n_controls + 2;
        idx_uz = N*n_states + (k-1)*n_controls + 3;
        
        model.lb(idx_ux) = -mission.u_max_x;
        model.ub(idx_ux) = mission.u_max_x;
        
        model.lb(idx_uy) = -mission.u_max_y;
        model.ub(idx_uy) = mission.u_max_y;
        
        model.lb(idx_uz) = -mission.u_max_z;
        model.ub(idx_uz) = mission.u_max_z;
    end
    
    % 등식 제약
    [Aeq, beq, ~, ~] = setupConstraints_QP(stm, N, RelInitState, RelFinalState, mission);
    model.A = Aeq;
    model.rhs = beq;
    model.sense = repmat('=', length(beq), 1);
    
    % Gurobi 파라미터
    params = struct();
    params.OutputFlag = 0;
    params.Method = 2;
    params.FeasibilityTol = 1e-6;
    params.OptimalityTol = 1e-6;
    
    % 해결
    result = gurobi(model, params);
    
    if strcmp(result.status, 'OPTIMAL')
        exitflag = 1;
        sol = result.x;
        
        X_opt = reshape(sol(1:N*n_states), n_states, N);
        U_opt = reshape(sol(N*n_states+1:end), n_controls, N-1);
    else
        exitflag = 0;
        X_opt = []; U_opt = [];
    end
end

function [Aeq, beq] = setupEqualityConstraints_LP(stm, N, RelInitState, RelFinalState, n_vars)
    % LP용 등식 제약 (보조 변수 제외)
    n_states = 6;
    n_controls = 3;
    n_x = N * n_states;
    
    % 등식 제약: 초기조건(6) + 동역학(6*(N-1)) + 최종조건(6)
    n_eq = n_states + n_states*(N-1) + n_states;
    
    % 비영 요소 개수 정확히 계산하여 메모리 효율성 향상
    nnz_est = n_states + (N-1)*(n_states*2 + n_controls) + n_states;
    Aeq = spalloc(n_eq, n_vars, nnz_est);
    beq = zeros(n_eq, 1);
    
    row_idx = 0;
    
    % 초기 조건
    Aeq(1:n_states, 1:n_states) = speye(n_states);
    beq(1:n_states) = RelInitState;
    row_idx = n_states;
    
    % 동역학 제약
    for k = 1:N-1
        rows = row_idx + (1:n_states);
        
        col_xkp1 = k*n_states + (1:n_states);
        Aeq(rows, col_xkp1) = speye(n_states);
        
        col_xk = (k-1)*n_states + (1:n_states);
        Aeq(rows, col_xk) = -sparse(stm.Ak(:,:,k));
        
        col_uk = n_x + (k-1)*n_controls + (1:n_controls);
        Aeq(rows, col_uk) = -sparse(stm.Bk(:,1:3,k));
        
        row_idx = row_idx + n_states;
    end
    
    % 최종 상태 조건
    rows = row_idx + (1:n_states);
    col_xN = (N-1)*n_states + (1:n_states);
    Aeq(rows, col_xN) = speye(n_states);
    beq(rows) = RelFinalState;
end

function [Aeq, beq, A_ineq, b_ineq] = setupConstraints_QP(stm, N, RelInitState, RelFinalState, mission)
    % QP용 제약 설정 (박스 제약은 Gurobi bounds로 처리)
    n_states = 6;
    n_controls = 3;
    n_vars = N*n_states + (N-1)*n_controls;
    
    % 등식 제약
    n_eq = n_states + n_states*(N-1) + n_states;
    Aeq = sparse(n_eq, n_vars);
    beq = zeros(n_eq, 1);
    
    row_idx = 0;
    
    Aeq(1:n_states, 1:n_states) = speye(n_states);
    beq(1:n_states) = RelInitState;
    row_idx = n_states;
    
    for k = 1:N-1
        rows = row_idx + (1:n_states);
        col_xkp1 = k*n_states + (1:n_states);
        Aeq(rows, col_xkp1) = speye(n_states);
        col_xk = (k-1)*n_states + (1:n_states);
        Aeq(rows, col_xk) = -sparse(stm.Ak(:,:,k));
        col_uk = N*n_states + (k-1)*n_controls + (1:n_controls);
        Aeq(rows, col_uk) = -sparse(stm.Bk(:,1:3,k));
        row_idx = row_idx + n_states;
    end
    
    rows = row_idx + (1:n_states);
    col_xN = (N-1)*n_states + (1:n_states);
    Aeq(rows, col_xN) = speye(n_states);
    beq(rows) = RelFinalState;
    
    % 부등식 제약 (QP에서는 필요하지만 Gurobi 모델에서는 bounds로 처리)
    A_ineq = [];
    b_ineq = [];
end

function total_dV = computeTotalDeltaV(U_opt, dt)
    total_dV = 0;
    for k = 1:size(U_opt, 2)
        total_dV = total_dV + norm(U_opt(:,k)) * dt;
    end
end

function total_dV = computeTotalDeltaV_BoxConstraints(U_opt, dt)
    % 박스 제약에서의 총 ΔV 계산 (각 축별 독립적)
    total_dV = 0;
    for k = 1:size(U_opt, 2)
        % L1 norm: |u_x| + |u_y| + |u_z|
        total_dV = total_dV + sum(abs(U_opt(:,k))) * dt;
    end
end

function visualizeComparison_SOCP(results_all, time_vec, mission, chief)
    figure('Position', [50, 50, 1800, 1000]);
    
    n_cases = length(results_all);
    colors = jet(n_cases);
    
    % 3D 궤적 비교
    subplot(2,4,1);
    hold on;
    for i = 1:n_cases
        X = results_all{i}.X_opt;
        plot3(X(1,:)/1000, X(3,:)/1000, X(5,:)/1000, '-', ...
              'Color', colors(i,:), 'LineWidth', 2);
    end
    plot3(0, 0, 0, 'k*', 'MarkerSize', 15);
    grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('3D 궤적 (LVLH)');
    view(45, 30);
    
    % 상대 거리
    subplot(2,4,2);
    hold on;
    for i = 1:n_cases
        X = results_all{i}.X_opt;
        rel_dist = sqrt(X(1,:).^2 + X(3,:).^2 + X(5,:).^2);
        semilogy(time_vec/3600, rel_dist, '-', 'Color', colors(i,:), 'LineWidth', 2);
    end
    yline(mission.rf, 'r--', '목표', 'LineWidth', 2);
    grid on;
    xlabel('시간 [hours]'); ylabel('거리 [m]');
    title('접근 프로파일');
    
    % 추력 프로파일 (각 축별)
    subplot(2,4,3);
    hold on;
    for i = 1:n_cases
        U = results_all{i}.U_opt;
        plot(time_vec(1:end-1)/3600, U(1,:)*1000, 'r-', 'LineWidth', 1.5);
        plot(time_vec(1:end-1)/3600, U(2,:)*1000, 'g-', 'LineWidth', 1.5);
        plot(time_vec(1:end-1)/3600, U(3,:)*1000, 'b-', 'LineWidth', 1.5);
    end
    yline(mission.u_max_x*1000, 'r--', 'LineWidth', 1);
    yline(-mission.u_max_x*1000, 'r--', 'LineWidth', 1);
    grid on;
    xlabel('시간 [hours]'); ylabel('추력 [mm/s²]');
    title('추력 프로파일 (X:빨강, Y:초록, Z:파랑)');
    legend('u_x', 'u_y', 'u_z', 'Location', 'best');
    
    % 누적 ΔV
    subplot(2,4,4);
    hold on;
    dt = mean(diff(time_vec));
    for i = 1:n_cases
        U = results_all{i}.U_opt;
        dV_cum = zeros(1, size(U,2));
        for k = 2:size(U,2)
            % L1 norm for box constraints
            dV_cum(k) = dV_cum(k-1) + sum(abs(U(:,k-1)))*dt;
        end
        plot(time_vec(1:end-1)/3600, dV_cum, '-', 'Color', colors(i,:), 'LineWidth', 2);
    end
    grid on;
    xlabel('시간 [hours]'); ylabel('누적 ΔV [m/s]');
    title('연료 소비');
    
    % 추력 사용 히스토그램
    subplot(2,4,5);
    for i = 1:n_cases
        U = results_all{i}.U_opt;
        histogram(abs(U(1,:))*1000, 20, 'FaceColor', 'r', 'FaceAlpha', 0.3);
        hold on;
        histogram(abs(U(2,:))*1000, 20, 'FaceColor', 'g', 'FaceAlpha', 0.3);
        histogram(abs(U(3,:))*1000, 20, 'FaceColor', 'b', 'FaceAlpha', 0.3);
    end
    xlabel('추력 크기 [mm/s²]'); ylabel('빈도');
    title('추력 사용 분포');
    legend('X축', 'Y축', 'Z축');
    
    % 각 축별 총 ΔV
    subplot(2,4,6);
    for i = 1:n_cases
        U = results_all{i}.U_opt;
        dV_x = sum(abs(U(1,:))) * dt;
        dV_y = sum(abs(U(2,:))) * dt;
        dV_z = sum(abs(U(3,:))) * dt;
        
        bar([1 2 3], [dV_x, dV_y, dV_z]);
        set(gca, 'XTickLabel', {'X', 'Y', 'Z'});
        ylabel('ΔV [m/s]');
        title('축별 연료 소비');
    end
    
    % 최종 위치 확대
    subplot(2,4,7);
    hold on;
    for i = 1:n_cases
        X = results_all{i}.X_opt;
        plot(X(1,end), X(3,end), 'o', 'Color', colors(i,:), ...
             'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
    end
    plot(14.14, 14.14, 'k*', 'MarkerSize', 15);  % 목표 위치
    grid on;
    xlabel('X [m]'); ylabel('Y [m]');
    title('최종 위치 (확대)');
    axis equal;
    
    % 요약 테이블
    subplot(2,4,8);
    axis off;
    text(0.1, 0.9, '=== 결과 요약 ===', 'FontSize', 14, 'FontWeight', 'bold');
    
    text(0.1, 0.8, '방법', 'FontWeight', 'bold');
    text(0.4, 0.8, 'ΔV [m/s]', 'FontWeight', 'bold');
    text(0.7, 0.8, '오차 [m]', 'FontWeight', 'bold');
    
    for i = 1:min(n_cases, 6)
        y_pos = 0.7 - 0.1*(i-1);
        if isfield(results_all{i}, 'method')
            text(0.1, y_pos, results_all{i}.method);
        else
            text(0.1, y_pos, 'LP');
        end
        text(0.4, y_pos, sprintf('%.3f', results_all{i}.total_dV));
        text(0.7, y_pos, sprintf('%.2e', results_all{i}.final_error));
    end
    
    sgtitle('Fuel Optimal Rendezvous with Real Thruster Model (LP) - LVLH Frame');
end

%% 향상된 ECI 좌표계 변환 및 시각화 함수들

function [r_eci, v_eci] = keplerian2ECI_improved(elements, mu, t)
    % 개선된 케플러 궤도 요소를 ECI 좌표로 변환
    % Newton-Raphson 방법으로 더 정확한 케플러 방정식 해법
    
    a = elements(1);
    e = elements(2);
    i = elements(3);
    RAAN = elements(4);
    w = elements(5);
    M0 = elements(6);
    
    n = sqrt(mu/a^3);
    N = length(t);
    
    r_eci = zeros(3, N);
    v_eci = zeros(3, N);
    
    for k = 1:N
        % Mean anomaly at time t
        M = mod(M0 + n*t(k), 2*pi);
        
        % Solve Kepler's equation using Newton-Raphson method
        % Better initial guess for high eccentricity
        if e < 0.8
            E = M + e*sin(M);
        else
            E = pi;
        end
        
        % Newton-Raphson iteration
        tol = 1e-12;
        max_iter = 50;
        for iter = 1:max_iter
            f = E - e*sin(E) - M;
            fp = 1 - e*cos(E);
            E_new = E - f/fp;
            
            if abs(E_new - E) < tol
                break;
            end
            E = E_new;
        end
        
        % True anomaly
        cosE = cos(E);
        sinE = sin(E);
        sqrt1pe = sqrt(1 + e);
        sqrt1me = sqrt(1 - e);
        
        nu = 2*atan2(sqrt1pe*sinE, sqrt1me*(1 + cosE));
        
        % Orbit radius
        r = a*(1 - e*cosE);
        
        % Position and velocity in perifocal frame
        cosnu = cos(nu);
        sinnu = sin(nu);
        
        r_pqw = r * [cosnu; sinnu; 0];
        
        % Velocity in perifocal frame
        p = a*(1 - e^2);  % Semi-latus rectum
        h = sqrt(mu*p);   % Specific angular momentum
        
        v_pqw = [-mu/h * sinnu; 
                 mu/h * (e + cosnu); 
                 0];
        
        % Rotation matrices
        cosRAAN = cos(RAAN);
        sinRAAN = sin(RAAN);
        cosi = cos(i);
        sini = sin(i);
        cosw = cos(w);
        sinw = sin(w);
        
        % Combined rotation matrix (PQW to ECI)
        R11 = cosRAAN*cosw - sinRAAN*sinw*cosi;
        R12 = -cosRAAN*sinw - sinRAAN*cosw*cosi;
        R13 = sinRAAN*sini;
        R21 = sinRAAN*cosw + cosRAAN*sinw*cosi;
        R22 = -sinRAAN*sinw + cosRAAN*cosw*cosi;
        R23 = -cosRAAN*sini;
        R31 = sinw*sini;
        R32 = cosw*sini;
        R33 = cosi;
        
        R_pqw2eci = [R11, R12, R13;
                     R21, R22, R23;
                     R31, R32, R33];
        
        % Transform to ECI
        r_eci(:,k) = R_pqw2eci * r_pqw;
        v_eci(:,k) = R_pqw2eci * v_pqw;
    end
end

function plotEarth3D_textured(ax, Req)
    % 텍스처가 있는 3D 지구 그리기
    
    % 구체 생성 (더 높은 해상도)
    [X, Y, Z] = sphere(100);
    X = X * Req/1000;  % km 단위
    Y = Y * Req/1000;
    Z = Z * Req/1000;
    
    % 지구 텍스처 로드 시도
    try
        % MATLAB에 포함된 지구 텍스처 사용
        load('topo.mat', 'topo', 'topomap1');
        
        % 텍스처를 구체에 매핑
        earth_surface = surf(ax, X, Y, Z);
        set(earth_surface, 'FaceColor', 'texturemap', ...
            'CData', topo, 'EdgeColor', 'none', 'FaceAlpha', 1);
        colormap(ax, topomap1);
        
    catch
        % 텍스처를 로드할 수 없는 경우 기본 색상 사용
        fprintf('지구 텍스처를 로드할 수 없습니다. 기본 색상을 사용합니다.\n');
        
        % 기본 파란색 지구
        surf(ax, X, Y, Z, 'FaceColor', [0.2 0.4 0.8], ...
             'EdgeColor', 'none', 'FaceAlpha', 0.9);
        
        % 간단한 대륙 표현
        hold(ax, 'on');
        
        % 적도
        theta = linspace(0, 2*pi, 200);
        plot3(ax, Req/1000*cos(theta), Req/1000*sin(theta), zeros(size(theta)), ...
              'y-', 'LineWidth', 1.5);
        
        % 주요 경도선
        for lon = 0:60:300
            phi = linspace(-pi/2, pi/2, 100);
            x = Req/1000 * cos(phi) * cos(deg2rad(lon));
            y = Req/1000 * cos(phi) * sin(deg2rad(lon));
            z = Req/1000 * sin(phi);
            plot3(ax, x, y, z, 'k-', 'LineWidth', 0.3, 'Color', [0.6 0.6 0.6]);
        end
        
        % 주요 위도선
        for lat = [-66.5, -23.5, 23.5, 66.5]  % 회귀선과 극권
            theta = linspace(0, 2*pi, 200);
            x = Req/1000 * cos(deg2rad(lat)) * cos(theta);
            y = Req/1000 * cos(deg2rad(lat)) * sin(theta);
            z = Req/1000 * sin(deg2rad(lat)) * ones(size(theta));
            plot3(ax, x, y, z, 'k--', 'LineWidth', 0.5, 'Color', [0.7 0.7 0.7]);
        end
    end
    
    % 축 설정
    xlabel(ax, 'X [km]');
    ylabel(ax, 'Y [km]');
    zlabel(ax, 'Z [km]');
    
    % 조명 효과
    lighting gouraud;
    light('Position', [1 0 1], 'Style', 'infinite');
    light('Position', [-1 -1 -1], 'Style', 'infinite', 'Color', [0.3 0.3 0.3]);
    
    axis(ax, 'equal');
    grid(ax, 'on');
    grid(ax, 'minor');
    
    % 배경색 설정 (우주)
    set(ax, 'Color', [0.05 0.05 0.1]);
    set(ax, 'GridColor', [0.3 0.3 0.3]);
    set(ax, 'MinorGridColor', [0.2 0.2 0.2]);
end

function visualizeECITrajectories_enhanced(results_all, time_vec, mission, chief, params)
    % 향상된 ECI 좌표계 시각화 (수정된 버전)
    
    figure('Position', [50, 50, 1800, 900]);
    
    % Chief 위성의 ECI 궤적 계산 (개선된 함수 사용)
    [r_chief_eci, v_chief_eci] = keplerian2ECI_improved(chief.elements, params.mu, time_vec);
    
    n_cases = length(results_all);
    colors = jet(n_cases);
    
    %% 서브플롯 1: 전체 궤도 보기 (3D)
    subplot(2, 3, [1, 4]);
    ax1 = gca;
    
    % 텍스처가 있는 지구 그리기
    plotEarth3D_textured(ax1, params.Req);
    hold on;
    
    % Chief 궤도 (전체)
    plot3(ax1, r_chief_eci(1,:)/1000, r_chief_eci(2,:)/1000, r_chief_eci(3,:)/1000, ...
          'w-', 'LineWidth', 2.5, 'DisplayName', 'Chief');
    
    % Deputy 궤적들
    for i = 1:n_cases
        X_rel = results_all{i}.X_opt;
        
        % 상대 위치를 ECI로 변환
        r_deputy_eci = zeros(3, size(X_rel, 2));
        v_deputy_eci = zeros(3, size(X_rel, 2));
        
        for k = 1:size(X_rel, 2)
            r_rel_lvlh = [X_rel(1,k); X_rel(3,k); X_rel(5,k)];
            v_rel_lvlh = [X_rel(2,k); X_rel(4,k); X_rel(6,k)];
            
            % 수정된 변환 사용
            [r_rel_eci, v_rel_eci] = lvlh2eci_improved(r_rel_lvlh, v_rel_lvlh, ...
                                                       r_chief_eci(:,k), v_chief_eci(:,k));
            r_deputy_eci(:,k) = r_chief_eci(:,k) + r_rel_eci;
            v_deputy_eci(:,k) = v_chief_eci(:,k) + v_rel_eci;
        end
        
        if isfield(results_all{i}, 'lambda')
            label = sprintf('Deputy λ=%.0e', results_all{i}.lambda);
        else
            label = sprintf('Deputy %s', results_all{i}.method);
        end
        
        plot3(ax1, r_deputy_eci(1,:)/1000, r_deputy_eci(2,:)/1000, r_deputy_eci(3,:)/1000, ...
              '-', 'Color', colors(i,:), 'LineWidth', 2, 'DisplayName', label);
    end
    
    % 시작점과 끝점 표시
    plot3(ax1, r_chief_eci(1,1)/1000, r_chief_eci(2,1)/1000, r_chief_eci(3,1)/1000, ...
          'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
    plot3(ax1, r_chief_eci(1,end)/1000, r_chief_eci(2,end)/1000, r_chief_eci(3,end)/1000, ...
          'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'End');
    
    title(ax1, 'ECI Frame - Full Orbit View', 'Color', 'w');
    legend(ax1, 'Location', 'northeastoutside', 'TextColor', 'w', 'Color', [0.1 0.1 0.2]);
    view(ax1, 45, 30);
    set(ax1, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
    
    %% 서브플롯 2: XY 평면 투영
    subplot(2, 3, 2);
    ax2 = gca;
    hold on;
    
    % 지구 원
    theta = linspace(0, 2*pi, 100);
    plot(ax2, params.Req/1000*cos(theta), params.Req/1000*sin(theta), 'b-', 'LineWidth', 2);
    fill(ax2, params.Req/1000*cos(theta), params.Req/1000*sin(theta), [0.2 0.4 0.8], 'FaceAlpha', 0.3);
    
    % 궤도 투영
    plot(ax2, r_chief_eci(1,:)/1000, r_chief_eci(2,:)/1000, 'k-', 'LineWidth', 2);
    
    for i = 1:n_cases
        X_rel = results_all{i}.X_opt;
        r_deputy_eci = zeros(3, size(X_rel, 2));
        for k = 1:size(X_rel, 2)
            r_rel_lvlh = [X_rel(1,k); X_rel(3,k); X_rel(5,k)];
            v_rel_lvlh = [X_rel(2,k); X_rel(4,k); X_rel(6,k)];
            [r_rel_eci, ~] = lvlh2eci_improved(r_rel_lvlh, v_rel_lvlh, ...
                                               r_chief_eci(:,k), v_chief_eci(:,k));
            r_deputy_eci(:,k) = r_chief_eci(:,k) + r_rel_eci;
        end
        plot(ax2, r_deputy_eci(1,:)/1000, r_deputy_eci(2,:)/1000, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    
    xlabel(ax2, 'X [km]'); ylabel(ax2, 'Y [km]');
    title(ax2, 'XY Plane Projection');
    axis(ax2, 'equal');
    grid(ax2, 'on');
    
    %% 서브플롯 3: XZ 평면 투영
    subplot(2, 3, 3);
    ax3 = gca;
    hold on;
    
    % 지구 원
    plot(ax3, params.Req/1000*cos(theta), params.Req/1000*sin(theta), 'b-', 'LineWidth', 2);
    fill(ax3, params.Req/1000*cos(theta), params.Req/1000*sin(theta), [0.2 0.4 0.8], 'FaceAlpha', 0.3);
    
    % 궤도 투영
    plot(ax3, r_chief_eci(1,:)/1000, r_chief_eci(3,:)/1000, 'k-', 'LineWidth', 2);
    
    for i = 1:n_cases
        X_rel = results_all{i}.X_opt;
        r_deputy_eci = zeros(3, size(X_rel, 2));
        for k = 1:size(X_rel, 2)
            r_rel_lvlh = [X_rel(1,k); X_rel(3,k); X_rel(5,k)];
            v_rel_lvlh = [X_rel(2,k); X_rel(4,k); X_rel(6,k)];
            [r_rel_eci, ~] = lvlh2eci_improved(r_rel_lvlh, v_rel_lvlh, ...
                                               r_chief_eci(:,k), v_chief_eci(:,k));
            r_deputy_eci(:,k) = r_chief_eci(:,k) + r_rel_eci;
        end
        plot(ax3, r_deputy_eci(1,:)/1000, r_deputy_eci(3,:)/1000, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    end
    
    xlabel(ax3, 'X [km]'); ylabel(ax3, 'Z [km]');
    title(ax3, 'XZ Plane Projection');
    axis(ax3, 'equal');
    grid(ax3, 'on');
    
    %% 서브플롯 4: 접근 구간 확대 (3D)
    subplot(2, 3, [5, 6]);
    ax4 = gca;
    
    % 마지막 1/4 구간만 표시
    idx_start = round(3*length(time_vec)/4);
    
    hold on;
    set(ax4, 'Color', [0.05 0.05 0.1]);
    
    % Chief 궤도 (접근 구간)
    plot3(ax4, r_chief_eci(1,idx_start:end)/1000, ...
          r_chief_eci(2,idx_start:end)/1000, ...
          r_chief_eci(3,idx_start:end)/1000, ...
          'w-', 'LineWidth', 3, 'DisplayName', 'Chief');
    
    % Deputy 궤적들 (접근 구간)
    for i = 1:n_cases
        X_rel = results_all{i}.X_opt;
        
        r_deputy_eci = zeros(3, size(X_rel, 2));
        for k = 1:size(X_rel, 2)
            r_rel_lvlh = [X_rel(1,k); X_rel(3,k); X_rel(5,k)];
            v_rel_lvlh = [X_rel(2,k); X_rel(4,k); X_rel(6,k)];
            [r_rel_eci, ~] = lvlh2eci_improved(r_rel_lvlh, v_rel_lvlh, ...
                                               r_chief_eci(:,k), v_chief_eci(:,k));
            r_deputy_eci(:,k) = r_chief_eci(:,k) + r_rel_eci;
        end
        
        if isfield(results_all{i}, 'lambda')
            label = sprintf('Deputy λ=%.0e', results_all{i}.lambda);
        else
            label = sprintf('Deputy %s', results_all{i}.method);
        end
        
        plot3(ax4, r_deputy_eci(1,idx_start:end)/1000, ...
              r_deputy_eci(2,idx_start:end)/1000, ...
              r_deputy_eci(3,idx_start:end)/1000, ...
              '-', 'Color', colors(i,:), 'LineWidth', 2.5, ...
              'DisplayName', label);
    end
    
    % 최종 위치 강조
    plot3(ax4, r_chief_eci(1,end)/1000, r_chief_eci(2,end)/1000, r_chief_eci(3,end)/1000, ...
          'r*', 'MarkerSize', 20, 'DisplayName', 'Target');
    
    % 축 설정
    xlabel(ax4, 'X [km]'); ylabel(ax4, 'Y [km]'); zlabel(ax4, 'Z [km]');
    title(ax4, 'Approach Phase - Zoomed View', 'Color', 'w');
    legend(ax4, 'Location', 'best', 'TextColor', 'w', 'Color', [0.1 0.1 0.2]);
    grid(ax4, 'on');
    grid(ax4, 'minor');
    axis(ax4, 'equal');
    view(ax4, 45, 30);
    set(ax4, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
    set(ax4, 'GridColor', [0.3 0.3 0.3]);
    set(ax4, 'MinorGridColor', [0.2 0.2 0.2]);
    
    % 전체 그림 제목
    sgtitle('Enhanced ECI Frame Trajectory Visualization', 'Color', 'w');
    set(gcf, 'Color', [0.1 0.1 0.15]);
end

function [r_rel_eci, v_rel_eci] = lvlh2eci_improved(r_lvlh, v_lvlh, r_chief_eci, v_chief_eci)
    % LVLH to ECI 변환 (위치와 속도 모두)
    
    r = norm(r_chief_eci);
    h = cross(r_chief_eci, v_chief_eci);
    
    % LVLH 프레임의 단위 벡터들 (ECI 기준)
    R_hat = r_chief_eci / r;
    W_hat = h / norm(h);
    S_hat = cross(W_hat, R_hat);
    
    % LVLH to ECI 회전 행렬
    R_lvlh2eci = [R_hat, S_hat, W_hat];
    
    % 1. 상대 위치 변환 (ECI 기준)
    r_rel_eci = R_lvlh2eci * r_lvlh;
    
    % 2. LVLH 프레임의 각속도 벡터 (ECI 기준)
    omega_eci = h / r^2;
    
    % 3. ECI 기준 상대 속도 계산 (수송 정리 적용)
    % v_rel_eci = (ECI로 회전된 v_rel_lvlh) + (프레임 회전으로 인한 속도)
    v_rel_eci = (R_lvlh2eci * v_lvlh) + cross(omega_eci, r_rel_eci);
end