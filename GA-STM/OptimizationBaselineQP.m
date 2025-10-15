%% 연료 최적 랑데부 - Hard Constraint 버전
% 최종 상태를 hard constraint로 설정하여 정확한 도달 보장
clear; clc; close all;

%% 1. 시스템 파라미터 설정
params = struct();
params.Req = 6378.137e3;           % 지구 반지름 [m]
params.mu = 3.986004415e14;        % 중력 상수 [m^3/s^2]
params.J2 = 1082.63e-6;            % J2 섭동 계수
params.tol = 1e-12;                % 케플러 방정식 허용 오차
params.safetyAltitude = 100e3;     % 안전 고도 [m]

% 주위성 궤도 요소
chief = struct();
chief.altitude = 700e3;            % 고도 [m]
chief.a = params.Req + chief.altitude;
chief.e = 0.001;
chief.i = deg2rad(98.2);
chief.RAAN = deg2rad(30);
chief.w = deg2rad(0);
chief.M = deg2rad(60);

chief.elements = [chief.a, chief.e, chief.i, chief.RAAN, chief.w, chief.M];
chief.n = sqrt(params.mu/chief.a^3);
chief.T = 2*pi/chief.n;

%% 2. 미션 파라미터
mission = struct();
mission.r0 = 50000;                % 초기 거리 [m]
mission.rf = 20;                   % 최종 거리 [m]
mission.u_max = 0.1;               % 최대 추력 [m/s^2]
mission.constraint_type = 'hard';   % 'hard' or 'soft'
mission.final_state_tol = 1e-3;    % Hard constraint 허용 오차 [m]

% 초기/최종 상태: custom_states.csv가 있으면 이를 우선 사용
% 다양한 상대경로에서 찾도록 후보 경로를 구성
here = fileparts(mfilename('fullpath'));
csv_candidates = { ...
    fullfile(pwd,'custom_states.csv'), ...
    fullfile(here,'custom_states.csv'), ...
    fullfile(here,'..','custom_states.csv') };
csv_path = '';
for ii=1:numel(csv_candidates)
    if exist(csv_candidates{ii},'file')
        csv_path = csv_candidates{ii};
        break;
    end
end
if ~isempty(csv_path)
    Xio = readmatrix(csv_path);
    RelInitState = Xio(1,:).';
    RelFinalState = Xio(2,:).';
    fprintf('Using custom start/end from: %s\n', csv_path);
    fprintf('  Start [x vx y vy z vz] (m): % .3f % .3f % .3f % .3f % .3f % .3f\n', RelInitState);
    fprintf('  Goal  [x vx y vy z vz] (m): % .3f % .3f % .3f % .3f % .3f % .3f\n', RelFinalState);
else
    [RelInitState, RelFinalState] = setupInitialFinalStates(mission.r0, mission.rf);
    fprintf('Using default start/end from setupInitialFinalStates().\n');
end

% 시간 설정
time_opts = struct();
time_opts.total_orbits = 1;
time_opts.dt = 60;
total_time = time_opts.total_orbits * chief.T;
time_vec = 0:time_opts.dt:total_time;
N = length(time_vec);

fprintf('=== 연료 최적 랑데부 (Hard Constraint) ===\n');
init_range = norm(RelInitState([1,3,5]));
final_range = norm(RelFinalState([1,3,5]));
fprintf('초기 거리(입력 기반): %.3f km\n', init_range/1000);
fprintf('최종 거리(입력 기반): %.3f m (Hard Constraint)\n', final_range);
fprintf('최대 추력: %.1f mm/s²\n', mission.u_max*1000);
fprintf('시간 스텝: %d\n\n', N);

%% 3. STM 계산
fprintf('STM 계산 중...\n');
tic;
stm = computeSTM(params, chief, time_vec, RelInitState);
time_stm = toc;
fprintf('  완료: %.3f 초\n\n', time_stm);

%% 4. 최적화 문제 설정 - HARD CONSTRAINT
fprintf('최적화 문제 구성 중 (Hard Constraint)...\n');

% 다양한 연료 가중치 테스트
R_values = [1e4]; % fixed R
results_all = cell(length(R_values), 1);

for r_idx = 1:length(R_values)
    weights.R = R_values(r_idx);
    
    fprintf('\n--- R = %.0e ---\n', weights.R);
    
    % QP 문제 구성
    tic;
    [H, f, Aeq, beq, A_ineq, b_ineq] = setupOptimizationProblemHard(stm, N, ...
        RelInitState, RelFinalState, mission, weights);
    time_setup = toc;
    
    % 문제 크기
    n_vars = size(H, 1);
    fprintf('문제 구성 완료: %.3f 초\n', time_setup);
    fprintf('  결정 변수: %d\n', n_vars);
    fprintf('  등식 제약: %d (최종 상태 포함)\n', size(Aeq, 1));
    fprintf('  부등식 제약: %d\n', size(A_ineq, 1));
    
    % 초기 추정치
    x0 = generateInitialGuess(stm, N, RelInitState, RelFinalState, mission);
    
    % 최적화 실행
    fprintf('최적화 실행 중...\n');
    tic;
    
    if exist('quadprog', 'file')
        options = optimoptions('quadprog', ...
            'Display', 'off', ...
            'MaxIterations', 3000, ...
            'OptimalityTolerance', 1e-10, ...
            'ConstraintTolerance', 1e-10);
        
        [sol, fval, exitflag, output] = quadprog(H, f, A_ineq, b_ineq, ...
                                                 Aeq, beq, [], [], x0, options);
    else
        error('quadprog not available');
    end
    
    time_solve = toc;
    
    % 결과 추출
    [X_opt, U_opt] = extractSolution(sol, N);
    total_dV = computeTotalDeltaV(U_opt, time_opts.dt);
    
    % 최종 상태 정확도 확인
    final_error = norm(X_opt(:,end) - RelFinalState);
    final_pos_error = norm(X_opt([1,3,5],end) - RelFinalState([1,3,5]));
    final_vel_error = norm(X_opt([2,4,6],end) - RelFinalState([2,4,6]));
    
    fprintf('  최적화 완료: %.3f 초 (exitflag: %d)\n', time_solve, exitflag);
    fprintf('  총 ΔV: %.3f m/s\n', total_dV);
    fprintf('  최종 위치 오차: %.2e m\n', final_pos_error);
    fprintf('  최종 속도 오차: %.2e m/s\n', final_vel_error);
    
    % 결과 저장
    results_all{r_idx} = struct('R', weights.R, 'X_opt', X_opt, 'U_opt', U_opt, ...
                               'total_dV', total_dV, 'final_error', final_error, ...
                               'exitflag', exitflag);
end

%% 5. 결과 비교 시각화
visualizeComparison(results_all, time_vec, mission, chief);
% Persist latest QP result for downstream comparison
try
    save(fullfile(fileparts(mfilename('fullpath')),'qp_last.mat'), ...
         'X_opt','U_opt','time_vec','RelInitState','RelFinalState','weights');
catch
end

%% 6. Soft vs Hard Constraint 비교 (선택사항)
if false  % true로 변경하여 실행
    fprintf('\n\n=== Soft Constraint 비교 ===\n');
    weights_soft.R = 1;
    weights_soft.Q_pos = 1e6;
    weights_soft.Q_vel = 1e4;
    
    [H_soft, f_soft, Aeq_soft, beq_soft, A_ineq_soft, b_ineq_soft] = ...
        setupOptimizationProblemSoft(stm, N, RelInitState, RelFinalState, mission, weights_soft);
    
    [sol_soft, ~, ~, ~] = quadprog(H_soft, f_soft, A_ineq_soft, b_ineq_soft, ...
                                   Aeq_soft, beq_soft, [], [], x0, options);
    
    [X_soft, U_soft] = extractSolution(sol_soft, N);
    total_dV_soft = computeTotalDeltaV(U_soft, time_opts.dt);
    final_error_soft = norm(X_soft(:,end) - RelFinalState);
    
    fprintf('Soft Constraint 결과:\n');
    fprintf('  총 ΔV: %.3f m/s\n', total_dV_soft);
    fprintf('  최종 상태 오차: %.2e m\n', final_error_soft);
end

%% 보조 함수들

function [RelInitState, RelFinalState] = setupInitialFinalStates(r0, rf)
    dx0 = r0 * 0.7071;
    dy0 = r0 * 0.7071;
    dz0 = 10000;
    RelInitState = [dx0; 0; dy0; 0; dz0; 0];
    
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

function [H, f, Aeq, beq, A_ineq, b_ineq] = setupOptimizationProblemHard(stm, N, ...
    RelInitState, RelFinalState, mission, weights)
    
    n_states = 6;
    n_controls = 3;
    n_vars = N*n_states + (N-1)*n_controls;
    
    % 목적 함수: 연료만 최소화 (최종 상태는 hard constraint로)
    H = sparse(n_vars, n_vars);
    f = sparse(n_vars, 1);
    
    % 연료 비용
    R = weights.R * eye(3);
    for k = 1:N-1
        idx = N*n_states + (k-1)*n_controls + (1:n_controls);
        H(idx, idx) = 2*R;
    end
    
    % Regularization
    H = H + 1e-10 * speye(n_vars);
    
    % 등식 제약
    % 크기: 초기조건(6) + 동역학(6*(N-1)) + 최종조건(6) = 6*N + 6
    n_eq = n_states + n_states*(N-1) + n_states;
    Aeq = sparse(n_eq, n_vars);
    beq = zeros(n_eq, 1);
    
    row_idx = 0;
    
    % 1. 초기 조건
    Aeq(1:n_states, 1:n_states) = speye(n_states);
    beq(1:n_states) = RelInitState;
    row_idx = n_states;
    
    % 2. 동역학 제약
    for k = 1:N-1
        rows = row_idx + (1:n_states);
        
        % x(k+1)
        col_xkp1 = k*n_states + (1:n_states);
        Aeq(rows, col_xkp1) = speye(n_states);
        
        % -A(k)*x(k)
        col_xk = (k-1)*n_states + (1:n_states);
        Aeq(rows, col_xk) = -sparse(stm.Ak(:,:,k));
        
        % -B(k)*u(k)
        col_uk = N*n_states + (k-1)*n_controls + (1:n_controls);
        Aeq(rows, col_uk) = -sparse(stm.Bk(:,1:3,k));
        
        row_idx = row_idx + n_states;
    end
    
    % 3. 최종 상태 조건 (HARD CONSTRAINT)
    rows = row_idx + (1:n_states);
    col_xN = (N-1)*n_states + (1:n_states);
    Aeq(rows, col_xN) = speye(n_states);
    beq(rows) = RelFinalState;
    
    % 부등식 제약 (추력 제한)
    n_ineq = 2*(N-1)*n_controls;
    A_ineq = sparse(n_ineq, n_vars);
    b_ineq = zeros(n_ineq, 1);
    
    row_idx = 0;
    for k = 1:N-1
        for j = 1:n_controls
            idx_u = N*n_states + (k-1)*n_controls + j;
            
            % u <= u_max
            row_idx = row_idx + 1;
            A_ineq(row_idx, idx_u) = 1;
            b_ineq(row_idx) = mission.u_max;
            
            % -u <= u_max
            row_idx = row_idx + 1;
            A_ineq(row_idx, idx_u) = -1;
            b_ineq(row_idx) = mission.u_max;
        end
    end
end

function [H, f, Aeq, beq, A_ineq, b_ineq] = setupOptimizationProblemSoft(stm, N, ...
    RelInitState, RelFinalState, mission, weights)
    % Soft constraint 버전 (비교용)
    
    n_states = 6;
    n_controls = 3;
    n_vars = N*n_states + (N-1)*n_controls;
    
    H = sparse(n_vars, n_vars);
    f = sparse(n_vars, 1);
    
    % 연료 비용
    R = weights.R * eye(3);
    for k = 1:N-1
        idx = N*n_states + (k-1)*n_controls + (1:n_controls);
        H(idx, idx) = 2*R;
    end
    
    % 최종 상태 비용 (SOFT)
    Q_final = zeros(6);
    Q_final(1:2:5, 1:2:5) = weights.Q_pos * eye(3);
    Q_final(2:2:6, 2:2:6) = weights.Q_vel * eye(3);
    
    idx_final = (N-1)*n_states + (1:n_states);
    H(idx_final, idx_final) = H(idx_final, idx_final) + 2*Q_final;
    f(idx_final) = -2*Q_final * RelFinalState;
    
    H = H + 1e-10 * speye(n_vars);
    
    % 등식 제약 (초기 조건 + 동역학만)
    n_eq = n_states + n_states*(N-1);
    Aeq = sparse(n_eq, n_vars);
    beq = zeros(n_eq, 1);
    
    % 초기 조건
    Aeq(1:n_states, 1:n_states) = speye(n_states);
    beq(1:n_states) = RelInitState;
    
    % 동역학
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
    
    % 부등식 제약
    n_ineq = 2*(N-1)*n_controls;
    A_ineq = sparse(n_ineq, n_vars);
    b_ineq = zeros(n_ineq, 1);
    
    row_idx = 0;
    for k = 1:N-1
        for j = 1:n_controls
            idx_u = N*n_states + (k-1)*n_controls + j;
            
            row_idx = row_idx + 1;
            A_ineq(row_idx, idx_u) = 1;
            b_ineq(row_idx) = mission.u_max;
            
            row_idx = row_idx + 1;
            A_ineq(row_idx, idx_u) = -1;
            b_ineq(row_idx) = mission.u_max;
        end
    end
end

function x0 = generateInitialGuess(stm, N, RelInitState, RelFinalState, mission)
    n_states = 6;
    n_controls = 3;
    n_vars = N*n_states + (N-1)*n_controls;
    x0 = zeros(n_vars, 1);
    
    % 상태: 자연 운동 + 선형 보간
    X_natural = zeros(n_states, N);
    X_natural(:,1) = RelInitState;
    
    for k = 1:N-1
        X_natural(:,k+1) = stm.Ak(:,:,k) * X_natural(:,k);
    end
    
    % 선형 보간으로 조정
    alpha = linspace(0, 1, N);
    for k = 1:N
        x0((k-1)*n_states + (1:n_states)) = (1-alpha(k))*X_natural(:,k) + alpha(k)*RelFinalState;
    end
    
    % 제어: Bang-bang 근사
    for k = 1:round(N/4)
        idx = N*n_states + (k-1)*n_controls + (1:n_controls);
        x0(idx) = -mission.u_max * 0.5 * [1; 0; 0];
    end
    
    for k = round(3*N/4):N-1
        idx = N*n_states + (k-1)*n_controls + (1:n_controls);
        x0(idx) = mission.u_max * 0.5 * [1; 0; 0];
    end
end

function [X_opt, U_opt] = extractSolution(sol, N)
    n_states = 6;
    n_controls = 3;
    
    X_opt = zeros(n_states, N);
    U_opt = zeros(n_controls, N-1);
    
    for k = 1:N
        X_opt(:,k) = sol((k-1)*n_states + (1:n_states));
    end
    
    for k = 1:N-1
        U_opt(:,k) = sol(N*n_states + (k-1)*n_controls + (1:n_controls));
    end
end

function total_dV = computeTotalDeltaV(U_opt, dt)
    total_dV = 0;
    for k = 1:size(U_opt, 2)
        total_dV = total_dV + norm(U_opt(:,k)) * dt;
    end
end

function visualizeComparison(results_all, time_vec, mission, chief)
    fig = figure('Position', [50, 50, 1800, 1000]);
    
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
    axis equal; % ensure equal scaling for trajectory axes
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('3D 궤적 비교');
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
    legend(arrayfun(@(i) sprintf('R=%.0e', results_all{i}.R), 1:n_cases, 'UniformOutput', false));
    
    % 추력 크기
    subplot(2,4,3);
    hold on;
    for i = 1:n_cases
        U = results_all{i}.U_opt;
        thrust_mag = sqrt(sum(U.^2, 1));
        stairs(time_vec(1:end-1)/3600, thrust_mag*1000, '-', ...
               'Color', colors(i,:), 'LineWidth', 2);
    end
    yline(mission.u_max*1000, 'k--', '최대', 'LineWidth', 2);
    grid on;
    xlabel('시간 [hours]'); ylabel('추력 크기 [mm/s²]');
    title('추력 프로파일');
    
    % 누적 ΔV
    subplot(2,4,4);
    hold on;
    dt = mean(diff(time_vec));
    for i = 1:n_cases
        U = results_all{i}.U_opt;
        dV_cum = zeros(1, size(U,2));
        for k = 2:size(U,2)
            dV_cum(k) = dV_cum(k-1) + norm(U(:,k-1))*dt;
        end
        plot(time_vec(1:end-1)/3600, dV_cum, '-', 'Color', colors(i,:), 'LineWidth', 2);
    end
    grid on;
    xlabel('시간 [hours]'); ylabel('누적 ΔV [m/s]');
    title('연료 소비');
    
    % R vs 총 ΔV
    subplot(2,4,5);
    R_vals = arrayfun(@(i) results_all{i}.R, 1:n_cases);
    dV_vals = arrayfun(@(i) results_all{i}.total_dV, 1:n_cases);
    semilogx(R_vals, dV_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    grid on;
    xlabel('연료 가중치 R'); ylabel('총 ΔV [m/s]');
    title('R vs 연료 소비 (Hard Constraint)');
    
    % R vs 최종 오차
    subplot(2,4,6);
    error_vals = arrayfun(@(i) results_all{i}.final_error, 1:n_cases);
    loglog(R_vals, error_vals, 'ro-', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('연료 가중치 R'); ylabel('최종 상태 오차 [m]');
    title('R vs 최종 정확도');
    yline(mission.final_state_tol, 'g--', '허용 오차', 'LineWidth', 2);
    
    % 최종 위치 확대
    subplot(2,4,7);
    hold on;
    for i = 1:n_cases
        X = results_all{i}.X_opt;
        plot(X(1,end), X(3,end), 'o', 'Color', colors(i,:), ...
             'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
    end
    plot(results_all{1}.X_opt(1,end), results_all{1}.X_opt(3,end), 'k*', 'MarkerSize', 15);
    grid on;
    xlabel('X [m]'); ylabel('Y [m]');
    title('최종 위치 (확대)');
    axis equal;
    xlim([10, 20]); ylim([10, 20]);
    
    % 요약 테이블
    subplot(2,4,8);
    axis off;
    text(0.1, 0.9, '=== 결과 요약 ===', 'FontSize', 14, 'FontWeight', 'bold');
    text(0.1, 0.8, 'R 값', 'FontWeight', 'bold');
    text(0.4, 0.8, 'ΔV [m/s]', 'FontWeight', 'bold');
    text(0.7, 0.8, '오차 [m]', 'FontWeight', 'bold');
    
    for i = 1:min(n_cases, 6)
        y_pos = 0.7 - 0.1*(i-1);
        text(0.1, y_pos, sprintf('%.0e', results_all{i}.R));
        text(0.4, y_pos, sprintf('%.3f', results_all{i}.total_dV));
        text(0.7, y_pos, sprintf('%.2e', results_all{i}.final_error));
    end
    
    sgtitle('Hard Constraint 연료 최적화 - R 값에 따른 비교', 'FontSize', 16, 'FontWeight', 'bold');
    % Save figure when running in batch
    try
        print(fig, 'qp_comparison.png', '-dpng', '-r150');
    catch
    end
end
