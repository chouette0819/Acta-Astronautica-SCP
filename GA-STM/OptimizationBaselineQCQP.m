%% 에너지 최적 랑데부 - QCQP with Gurobi (L2 norm) + ECI 시각화 및 영상 저장
% 표준 QCQP 문제: L2 norm (에너지) 최소화
clear; clc; close all;

% Gurobi 경로 추가
addpath("C:\gurobi1202\win64\matlab")

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
mission.u_max = 0.05;              % 최대 추력 [m/s^2]

% 초기/최종 상태
[RelInitState, RelFinalState] = setupInitialFinalStates(mission.r0, mission.rf);

% 시간 설정
time_opts = struct();
time_opts.total_orbits = 3;        % 수치 안정성을 위해 3 orbit 사용
time_opts.dt = 60;                 % [s]
total_time = time_opts.total_orbits * chief.T;
time_vec = 0:time_opts.dt:total_time;
N = length(time_vec);

fprintf('=== 에너지 최적 랑데부 (QCQP L2) ===\n');
fprintf('초기 거리: %.1f km\n', mission.r0/1000);
fprintf('최종 거리: %.1f m\n', mission.rf);
fprintf('최대 추력: %.1f mm/s²\n', mission.u_max*1000);
fprintf('시간 스텝: %d\n', N);
fprintf('총 시간: %.2f hours\n\n', total_time/3600);

%% 3. STM 계산
fprintf('STM 계산 중...\n');
tic;
stm = computeSTM(params, chief, time_vec, RelInitState);
time_stm = toc;
fprintf('  완료: %.3f 초\n\n', time_stm);

%% 4. QCQP 문제 해결
fprintf('QCQP 문제 구성 및 해결 중...\n');

% 최적화 실행
[sol, info] = solveQCQP_L2(stm, N, RelInitState, RelFinalState, mission, time_opts);

fprintf('\n=== 최적화 결과 ===\n');
fprintf('상태: %s\n', info.status);
fprintf('해결 시간: %.3f 초\n', info.solve_time);
fprintf('에너지 (목적함수): %.6f\n', info.energy);
fprintf('총 ΔV: %.3f m/s\n', info.total_dV);
fprintf('최종 위치 오차: %.2e m\n', info.final_pos_error);
fprintf('최종 속도 오차: %.2e m/s\n', info.final_vel_error);

%% 5. 결과 시각화
visualizeResults(sol, info, time_vec, mission, chief, params);

%% 6. ECI 시각화 및 애니메이션 저장
fprintf('\n영상 생성 중...\n');

% 전체 뷰 애니메이션
fprintf('전체 뷰 영상 생성 중...\n');
createECIAnimation(sol, time_vec, chief, params, mission);

% Deputy 추적 뷰 애니메이션
fprintf('Deputy 추적 뷰 영상 생성 중...\n');
createDeputyFollowAnimation(sol, time_vec, chief, params, mission);

%% ==================== 함수 정의 ====================

function createDeputyFollowAnimation(sol, time_vec, chief, params, mission)
    % Deputy를 따라가는 카메라 뷰 애니메이션
    
    [r_chief_eci, r_deputy_eci] = computeECITrajectories(sol.X, time_vec, chief, params);
    
    % 비디오 설정
    video_filename = sprintf('rendezvous_deputy_follow_%.0fkm_to_%.0fm.mp4', ...
                           mission.r0/1000, mission.rf);
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;
    v.Quality = 95;
    open(v);
    
    % Figure 설정
    fig = figure('Position', [100, 100, 1280, 720], 'Color', 'k');
    ax = axes('Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
    hold on;
    
    % 지구 그리기 (출력 인수 없이)
    drawEarth(params.Req);  % earth_handle 제거
    
    % 궤도 경로 (전체, 매우 흐릿하게)
    plot3(r_chief_eci(1,:)/1e6, r_chief_eci(2,:)/1e6, r_chief_eci(3,:)/1e6, ...
          'Color', [0.2 0.2 0.3 0.2], 'LineWidth', 0.5);
    
    % 애니메이션 객체
    h_chief = plot3(0, 0, 0, 'd', 'Color', [0.7 0.7 1], ...
                    'MarkerSize', 8, 'MarkerFaceColor', [0.7 0.7 1]);
    h_deputy = plot3(0, 0, 0, 'o', 'Color', [1 0.5 0], ...
                     'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0]);
    h_line = plot3([0 0], [0 0], [0 0], 'y-', 'LineWidth', 2);
    
    % Deputy 궤적 (가까운 과거)
    h_trail = plot3(0, 0, 0, 'Color', [1 0.5 0 0.7], 'LineWidth', 3);
    
    % 추력 벡터 표시
    h_thrust = quiver3(0, 0, 0, 0, 0, 0, 'Color', [1 0 0], ...
                      'LineWidth', 2, 'MaxHeadSize', 0.5);
    
    % 정보 표시
    h_info = text(0.02, 0.98, '', 'Units', 'normalized', ...
                  'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top', ...
                  'BackgroundColor', [0 0 0 0.5], 'EdgeColor', 'w');
    
    % 축 설정
    axis equal;
    grid on;
    xlabel('X [×1000 km]', 'Color', 'w');
    ylabel('Y [×1000 km]', 'Color', 'w');
    zlabel('Z [×1000 km]', 'Color', 'w');
    title('Deputy 추적 카메라', 'Color', 'w', 'FontSize', 16);
    
    lighting gouraud;
    light('Position', [2 1 1], 'Color', [1 1 0.9]);
    light('Position', [-1 -1 -0.5], 'Color', [0.3 0.3 0.5]);
    
    % 별 배경 (고정)
    addStars(ax);
    
    % 애니메이션 프레임 생성
    skip = max(1, floor(length(time_vec)/300));
    trail_length = 50;
    
    try
        frame_count = 0;
        
        for k = 1:skip:length(time_vec)
            % 위성 위치 업데이트
            chief_pos = r_chief_eci(:,k)/1e6;
            deputy_pos = r_deputy_eci(:,k)/1e6;
            
            set(h_chief, 'XData', chief_pos(1), 'YData', chief_pos(2), 'ZData', chief_pos(3));
            set(h_deputy, 'XData', deputy_pos(1), 'YData', deputy_pos(2), 'ZData', deputy_pos(3));
            
            % 연결선
            set(h_line, 'XData', [chief_pos(1), deputy_pos(1)], ...
                        'YData', [chief_pos(2), deputy_pos(2)], ...
                        'ZData', [chief_pos(3), deputy_pos(3)]);
            
            % Deputy 궤적
            trail_start = max(1, k-trail_length*skip);
            trail_idx = trail_start:skip:k;
            if length(trail_idx) > 1
                set(h_trail, 'XData', r_deputy_eci(1,trail_idx)/1e6, ...
                             'YData', r_deputy_eci(2,trail_idx)/1e6, ...
                             'ZData', r_deputy_eci(3,trail_idx)/1e6);
            end
            
            % 추력 벡터 (크기 과장)
            if k < length(time_vec) && k <= size(sol.U, 2)
                thrust_scale = 2;  % 시각화를 위해 크게 확대
                u_vec = sol.U(:, k);
                
                % 현재 시간의 진근점이각 계산
                M = chief.M + chief.n * time_vec(k);
                E = solveKeplerEquation(M, chief.e);
                nu = 2 * atan2(sqrt(1 + chief.e) * sin(E/2), sqrt(1 - chief.e) * cos(E/2));
                
                % LVLH to ECI 변환
                [r_c, v_c] = oe2eci(chief.a, chief.e, chief.i, chief.RAAN, chief.w, nu, params.mu);
                R_lvlh2eci = getLVLH2ECIMatrix(r_c, v_c);
                u_eci = R_lvlh2eci * u_vec * thrust_scale;
                
                set(h_thrust, 'XData', deputy_pos(1), 'YData', deputy_pos(2), 'ZData', deputy_pos(3), ...
                             'UData', u_eci(1), 'VData', u_eci(2), 'WData', u_eci(3));
            end
            
            % 카메라 위치 설정 (Deputy 추적)
            % Deputy 속도 방향 계산
            if k < length(time_vec)
                deputy_vel = (r_deputy_eci(:,min(k+skip,end)) - r_deputy_eci(:,k)) / (skip*time_vec(2));
            else
                deputy_vel = (r_deputy_eci(:,k) - r_deputy_eci(:,max(1,k-skip))) / (skip*time_vec(2));
            end
            
            % 속도가 0인 경우 처리
            if norm(deputy_vel) < 1e-6
                deputy_vel = [1; 0; 0];  % 기본 방향
            end
            
            % 카메라 위치: Deputy 뒤쪽 상단
            cam_distance = 500;  % km
            cam_height = 200;    % km
            
            % 속도 방향의 반대
            vel_dir = deputy_vel / norm(deputy_vel);
            up_dir = deputy_pos / norm(deputy_pos);  % 지구 중심에서 바깥쪽
            right_dir = cross(vel_dir, up_dir);
            
            if norm(right_dir) < 1e-6  % 평행한 경우 처리
                right_dir = cross(vel_dir, [0; 0; 1]);
            end
            right_dir = right_dir / norm(right_dir);
            
            cam_pos = deputy_pos - vel_dir * cam_distance/1e3 + up_dir * cam_height/1e3;
            
            % 카메라 설정
            set(ax, 'CameraPosition', cam_pos);
            set(ax, 'CameraTarget', deputy_pos);
            set(ax, 'CameraUpVector', up_dir);
            set(ax, 'CameraViewAngle', 60);
            
            % 정보 업데이트
            rel_dist_k = norm(r_deputy_eci(:,k) - r_chief_eci(:,k));
            if k <= size(sol.U, 2)
                thrust_mag = norm(sol.U(:,k))*1000;
            else
                thrust_mag = 0;
            end
            
            info_str = sprintf(['시간: %.2f / %.2f hours\n' ...
                               '상대 거리: %.1f km\n' ...
                               '추력: %.3f mm/s²'], ...
                              time_vec(k)/3600, time_vec(end)/3600, ...
                              rel_dist_k/1000, thrust_mag);
            set(h_info, 'String', info_str);
            
            drawnow;
            
            % 프레임 캡처 및 저장
            frame = getframe(fig);
            if ~isempty(frame.cdata)
                writeVideo(v, frame);
                frame_count = frame_count + 1;
            end
        end
        
        % 마지막 프레임 유지 (2초)
        if frame_count > 0
            last_frame = getframe(fig);
            for i = 1:60
                writeVideo(v, last_frame);
            end
        end
        
    catch ME
        fprintf('애니메이션 생성 중 오류: %s\n', ME.message);
    end
    
    close(v);
    close(fig);
    
    if frame_count > 0
        fprintf('Deputy 추적 영상 저장 완료: %s (총 %d 프레임)\n', video_filename, frame_count);
    else
        fprintf('경고: 프레임이 생성되지 않았습니다.\n');
    end
end

function [RelInitState, RelFinalState] = setupInitialFinalStates(r0, rf)
    % 초기 및 최종 상태 설정
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
    % State Transition Matrix 계산
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

function [sol, info] = solveQCQP_L2(stm, N, RelInitState, RelFinalState, mission, time_opts)
    % L2 norm 최소화를 위한 QCQP 문제 해결
    
    tic;
    
    n_states = 6;
    n_controls = 3;
    dt = time_opts.dt;
    
    % 변수: [x(1)...x(N), u(1)...u(N-1)]
    n_x = N * n_states;
    n_u = (N-1) * n_controls;
    n_vars = n_x + n_u;
    
    % Gurobi 모델 생성
    model = struct();
    
    % 목적함수: minimize sum(||u_k||^2)
    model.Q = sparse(n_vars, n_vars);
    for k = 1:N-1
        idx_u = n_x + (k-1)*n_controls + (1:n_controls);
        for i = 1:n_controls
            model.Q(idx_u(i), idx_u(i)) = 1;
        end
    end
    model.obj = zeros(n_vars, 1);
    
    % 변수 bounds
    model.lb = -inf(n_vars, 1);
    model.ub = inf(n_vars, 1);
    
    % 선형 등식 제약: 초기조건 + 동역학 + 최종조건
    n_eq = n_states + n_states*(N-1) + n_states;
    A_eq = sparse(n_eq, n_vars);
    b_eq = zeros(n_eq, 1);
    
    row = 0;
    
    % 초기 조건
    A_eq(1:n_states, 1:n_states) = eye(n_states);
    b_eq(1:n_states) = RelInitState;
    row = n_states;
    
    % 동역학 제약
    for k = 1:N-1
        rows = row + (1:n_states);
        cols_xkp1 = k*n_states + (1:n_states);
        A_eq(rows, cols_xkp1) = eye(n_states);
        cols_xk = (k-1)*n_states + (1:n_states);
        A_eq(rows, cols_xk) = -stm.Ak(:,:,k);
        cols_uk = n_x + (k-1)*n_controls + (1:n_controls);
        A_eq(rows, cols_uk) = -stm.Bk(:,1:3,k);
        row = row + n_states;
    end
    
    % 최종 조건
    rows = row + (1:n_states);
    cols_xN = (N-1)*n_states + (1:n_states);
    A_eq(rows, cols_xN) = eye(n_states);
    b_eq(rows) = RelFinalState;
    
    model.A = A_eq;
    model.rhs = b_eq;
    model.sense = repmat('=', n_eq, 1);
    
    % Quadratic constraints: ||u(k)||_2^2 <= u_max^2
    model.quadcon = [];
    for k = 1:N-1
        qc = struct();
        qc.Qc = sparse(n_vars, n_vars);
        idx_u = n_x + (k-1)*n_controls + (1:n_controls);
        for i = 1:n_controls
            qc.Qc(idx_u(i), idx_u(i)) = 1;
        end
        qc.q = zeros(n_vars, 1);
        qc.rhs = mission.u_max^2;
        qc.sense = '<';
        model.quadcon = [model.quadcon, qc];
    end
    
    % Gurobi 파라미터 설정
    params_grb = struct();
    params_grb.OutputFlag = 0;
    params_grb.Method = 2;
    params_grb.BarConvTol = 1e-6;
    params_grb.FeasibilityTol = 1e-6;
    params_grb.OptimalityTol = 1e-6;
    params_grb.NumericFocus = 3;
    params_grb.ScaleFlag = 2;
    
    % 최적화 실행
    result = gurobi(model, params_grb);
    
    solve_time = toc;
    
    % 결과 추출
    if strcmp(result.status, 'OPTIMAL')
        x_opt = result.x;
        
        X_opt = zeros(n_states, N);
        U_opt = zeros(n_controls, N-1);
        
        for k = 1:N
            X_opt(:,k) = x_opt((k-1)*n_states + (1:n_states));
        end
        
        for k = 1:N-1
            U_opt(:,k) = x_opt(n_x + (k-1)*n_controls + (1:n_controls));
        end
        
        energy = 0;
        for k = 1:N-1
            energy = energy + U_opt(1,k)^2 + U_opt(2,k)^2 + U_opt(3,k)^2;
        end
        
        total_dV = 0;
        for k = 1:N-1
            total_dV = total_dV + norm(U_opt(:,k)) * dt;
        end
        
        final_error = norm(X_opt(:,end) - RelFinalState);
        final_pos_error = norm(X_opt([1,3,5],end) - RelFinalState([1,3,5]));
        final_vel_error = norm(X_opt([2,4,6],end) - RelFinalState([2,4,6]));
        
        sol = struct('X', X_opt, 'U', U_opt);
        info = struct('energy', energy, ...
                     'total_dV', total_dV, ...
                     'final_error', final_error, ...
                     'final_pos_error', final_pos_error, ...
                     'final_vel_error', final_vel_error, ...
                     'solve_time', solve_time, ...
                     'status', result.status);
    else
        error('최적화 실패: %s', result.status);
    end
end

function visualizeResults(sol, info, time_vec, mission, chief, params)
    % LVLH 결과 시각화
    figure('Position', [50, 50, 1400, 800]);
    
    X_opt = sol.X;
    U_opt = sol.U;
    
    % 1. 3D 궤적
    subplot(2,3,1);
    plot3(X_opt(1,:)/1000, X_opt(3,:)/1000, X_opt(5,:)/1000, 'b-', 'LineWidth', 2);
    hold on;
    plot3(0, 0, 0, 'k*', 'MarkerSize', 15);
    plot3(X_opt(1,1)/1000, X_opt(3,1)/1000, X_opt(5,1)/1000, 'go', ...
          'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot3(X_opt(1,end)/1000, X_opt(3,end)/1000, X_opt(5,end)/1000, 'ro', ...
          'MarkerSize', 10, 'MarkerFaceColor', 'r');
    grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    title('3D 상대 궤적 (LVLH)');
    view(45, 30);
    
    % 2. 상대 거리
    subplot(2,3,2);
    rel_dist = sqrt(X_opt(1,:).^2 + X_opt(3,:).^2 + X_opt(5,:).^2);
    semilogy(time_vec/3600, rel_dist, 'b-', 'LineWidth', 2);
    hold on;
    yline(mission.rf, 'r--', '목표', 'LineWidth', 2);
    grid on;
    xlabel('시간 [hours]'); ylabel('거리 [m]');
    title('상대 거리');
    
    % 3. 추력 크기
    subplot(2,3,3);
    thrust_mag = sqrt(sum(U_opt.^2, 1));
    stairs(time_vec(1:end-1)/3600, thrust_mag*1000, 'b-', 'LineWidth', 2);
    hold on;
    yline(mission.u_max*1000, 'k--', '최대', 'LineWidth', 2);
    grid on;
    xlabel('시간 [hours]'); ylabel('추력 크기 [mm/s²]');
    title('추력 프로파일');
    ylim([0, mission.u_max*1000*1.1]);
    
    % 4. 누적 ΔV
    subplot(2,3,4);
    dt = mean(diff(time_vec));
    dV_cum = zeros(1, size(U_opt,2));
    for k = 2:size(U_opt,2)
        dV_cum(k) = dV_cum(k-1) + norm(U_opt(:,k-1))*dt;
    end
    plot(time_vec(1:end-1)/3600, dV_cum, 'b-', 'LineWidth', 2);
    grid on;
    xlabel('시간 [hours]'); ylabel('누적 ΔV [m/s]');
    title(sprintf('연료 소비 (총 %.3f m/s)', info.total_dV));
    
    % 5. 추력 벡터 성분
    subplot(2,3,5);
    plot(time_vec(1:end-1)/3600, U_opt(1,:)*1000, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(time_vec(1:end-1)/3600, U_opt(2,:)*1000, 'g-', 'LineWidth', 1.5);
    plot(time_vec(1:end-1)/3600, U_opt(3,:)*1000, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('시간 [hours]'); ylabel('추력 성분 [mm/s²]');
    title('추력 벡터 성분');
    legend('u_x', 'u_y', 'u_z');
    
    % 6. 속도 성분
    subplot(2,3,6);
    plot(time_vec/3600, X_opt(2,:), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(time_vec/3600, X_opt(4,:), 'g-', 'LineWidth', 1.5);
    plot(time_vec/3600, X_opt(6,:), 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('시간 [hours]'); ylabel('속도 [m/s]');
    title('상대 속도 성분');
    legend('v_x', 'v_y', 'v_z');
    
    sgtitle(sprintf('에너지 최적 랑데부 (QCQP L2) - 에너지: %.3f, ΔV: %.3f m/s', ...
                    info.energy, info.total_dV), 'FontSize', 14);
                    
    % ECI 시각화 호출
    visualizeECI(sol, time_vec, chief, params);
end

function visualizeECI(sol, time_vec, chief, params)
    % ECI 좌표계에서 궤도 시각화
    
    [r_chief_eci, r_deputy_eci] = computeECITrajectories(sol.X, time_vec, chief, params);
    
    rel_distances = zeros(1, length(time_vec));
    for k = 1:length(time_vec)
        rel_distances(k) = norm(r_deputy_eci(:,k) - r_chief_eci(:,k));
    end
    
    figure('Position', [100, 100, 1400, 900], 'Color', 'k', 'Name', 'ECI 궤도 시각화');
    ax = axes('Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
    hold on;
    
    % 지구 그리기
    drawEarth(params.Req);
    
    % Chief 궤도
    plot3(r_chief_eci(1,:)/1e6, r_chief_eci(2,:)/1e6, r_chief_eci(3,:)/1e6, ...
          'Color', [0.4 0.6 1], 'LineWidth', 3, 'DisplayName', 'Chief (목표 위성)');
    
    % Deputy 궤도
    plot3(r_deputy_eci(1,:)/1e6, r_deputy_eci(2,:)/1e6, r_deputy_eci(3,:)/1e6, ...
          'Color', [1 0.6 0.2], 'LineWidth', 2.5, 'DisplayName', 'Deputy (추격 위성)');
    
    % 시작/종료 위치
    plot3(r_deputy_eci(1,1)/1e6, r_deputy_eci(2,1)/1e6, r_deputy_eci(3,1)/1e6, ...
          'o', 'Color', 'white', 'MarkerSize', 2, 'MarkerFaceColor', [0 1 0], ...
          'LineWidth', 3);
    plot3(r_deputy_eci(1,end)/1e6, r_deputy_eci(2,end)/1e6, r_deputy_eci(3,end)/1e6, ...
          'p', 'Color', 'white', 'MarkerSize', 2, 'MarkerFaceColor', [1 0 0], ...
          'LineWidth', 3);
    
    % 연결선
    n_lines = 10;
    for i = 1:n_lines
        idx = round(1 + (i-1)*(length(time_vec)-1)/(n_lines-1));
        plot3([r_chief_eci(1,idx), r_deputy_eci(1,idx)]/1e6, ...
              [r_chief_eci(2,idx), r_deputy_eci(2,idx)]/1e6, ...
              [r_chief_eci(3,idx), r_deputy_eci(3,idx)]/1e6, ...
              'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 1);
    end
    
    axis equal;
    grid on;
    xlabel('X [×1000 km]', 'Color', 'w');
    ylabel('Y [×1000 km]', 'Color', 'w');
    zlabel('Z [×1000 km]', 'Color', 'w');
    title(sprintf('우주 랑데부: %.1f km → %.1f m', ...
                  rel_distances(1)/1000, rel_distances(end)), ...
          'Color', 'w', 'FontSize', 16);
    
    lighting gouraud;
    light('Position', [2 1 1], 'Color', [1 1 0.9]);
    light('Position', [-1 -1 -0.5], 'Color', [0.3 0.3 0.5]);
    
    view(45, 20);
    legend('Location', 'best', 'TextColor', 'w', 'Color', 'none');
    
    addStars(ax);
end

function createECIAnimation(sol, time_vec, chief, params, mission)
    % ECI 좌표계 애니메이션 생성 및 저장
    
    [r_chief_eci, r_deputy_eci] = computeECITrajectories(sol.X, time_vec, chief, params);
    
    % 비디오 설정
    video_filename = sprintf('rendezvous_animation_%.0fkm_to_%.0fm.mp4', ...
                           mission.r0/1000, mission.rf);
    v = VideoWriter(video_filename, 'MPEG-4');
    v.FrameRate = 30;
    v.Quality = 95;
    open(v);
    
    % Figure 설정
    fig = figure('Position', [100, 100, 1280, 720], 'Color', 'k');
    ax = axes('Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
    hold on;
    
    % 지구 그리기
    drawEarth(params.Req);
    
    % 궤도 경로 (흐릿하게)
    plot3(r_chief_eci(1,:)/1e6, r_chief_eci(2,:)/1e6, r_chief_eci(3,:)/1e6, ...
          'Color', [0.3 0.3 0.5 0.3], 'LineWidth', 1);
    plot3(r_deputy_eci(1,:)/1e6, r_deputy_eci(2,:)/1e6, r_deputy_eci(3,:)/1e6, ...
          'Color', [0.5 0.3 0.1 0.3], 'LineWidth', 1);
    
    % 시작점 표시
    plot3(r_deputy_eci(1,1)/1e6, r_deputy_eci(2,1)/1e6, r_deputy_eci(3,1)/1e6, ...
          'o', 'Color', [0 1 0], 'MarkerSize', 8, 'MarkerFaceColor', [0 1 0]);
    
    % 애니메이션 객체
    h_chief = plot3(0, 0, 0, 'o', 'Color', [0.7 0.7 1], ...
                    'MarkerSize', 10, 'MarkerFaceColor', [0.7 0.7 1]);
    h_deputy = plot3(0, 0, 0, 'o', 'Color', [1 0.5 0], ...
                     'MarkerSize', 10, 'MarkerFaceColor', [1 0.5 0]);
    h_line = plot3([0 0], [0 0], [0 0], 'w-', 'LineWidth', 1.5);
    h_trail = plot3(0, 0, 0, 'Color', [1 0.5 0 0.5], 'LineWidth', 2);
    
    % 텍스트 정보
    h_time = text(0.02, 0.98, '', 'Units', 'normalized', ...
                  'Color', 'w', 'FontSize', 14, 'VerticalAlignment', 'top');
    h_dist = text(0.02, 0.92, '', 'Units', 'normalized', ...
                  'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
    
    % 축 설정
    axis equal;
    grid on;
    xlabel('X [×1000 km]', 'Color', 'w');
    ylabel('Y [×1000 km]', 'Color', 'w');
    zlabel('Z [×1000 km]', 'Color', 'w');
    title('우주 랑데부 애니메이션', 'Color', 'w', 'FontSize', 16);
    
    lighting gouraud;
    light('Position', [2 1 1], 'Color', [1 1 0.9]);
    
    % 축 범위 설정
    max_range = max([max(abs(r_chief_eci(:))), max(abs(r_deputy_eci(:)))]) / 1e6 * 1.2;
    xlim([-max_range, max_range]);
    ylim([-max_range, max_range]);
    zlim([-max_range, max_range]);
    
    addStars(ax);
    
    % 애니메이션 프레임 생성
    skip = max(1, floor(length(time_vec)/300));  % 약 300 프레임
    trail_length = 30;  % 궤적 길이
    
    for k = 1:skip:length(time_vec)
        % 위성 위치 업데이트
        set(h_chief, 'XData', r_chief_eci(1,k)/1e6, ...
                     'YData', r_chief_eci(2,k)/1e6, ...
                     'ZData', r_chief_eci(3,k)/1e6);
        set(h_deputy, 'XData', r_deputy_eci(1,k)/1e6, ...
                      'YData', r_deputy_eci(2,k)/1e6, ...
                      'ZData', r_deputy_eci(3,k)/1e6);
        
        % 연결선 업데이트
        set(h_line, 'XData', [r_chief_eci(1,k), r_deputy_eci(1,k)]/1e6, ...
                    'YData', [r_chief_eci(2,k), r_deputy_eci(2,k)]/1e6, ...
                    'ZData', [r_chief_eci(3,k), r_deputy_eci(3,k)]/1e6);
        
        % Deputy 궤적 업데이트
        trail_start = max(1, k-trail_length*skip);
        trail_idx = trail_start:skip:k;
        set(h_trail, 'XData', r_deputy_eci(1,trail_idx)/1e6, ...
                     'YData', r_deputy_eci(2,trail_idx)/1e6, ...
                     'ZData', r_deputy_eci(3,trail_idx)/1e6);
        
        % 텍스트 업데이트
        rel_dist_k = norm(r_deputy_eci(:,k) - r_chief_eci(:,k));
        set(h_time, 'String', sprintf('시간: %.2f / %.2f hours', ...
                                     time_vec(k)/3600, time_vec(end)/3600));
        set(h_dist, 'String', sprintf('상대 거리: %.1f km', rel_dist_k/1000));
        
        % 카메라 회전
        view(45 + k/length(time_vec)*180, 20);
        
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    % 마지막 프레임 추가 (2초)
    for i = 1:60
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    close(v);
    close(fig);
    
    fprintf('영상 저장 완료: %s\n', video_filename);
end

function drawEarth(R_earth)
    % NASA Blue Marble 텍스처를 사용한 실제 지구 그리기
    
    n_pts = 150;  % 더 높은 해상도
    [lon, lat] = meshgrid(linspace(-180, 180, n_pts), linspace(90, -90, n_pts));
    
    % 구면 좌표를 데카르트 좌표로 변환
    X = R_earth/1e6 * cosd(lat) .* cosd(lon);
    Y = R_earth/1e6 * cosd(lat) .* sind(lon);
    Z = R_earth/1e6 * sind(lat);
    
    try
        % NASA Blue Marble 이미지 다운로드 시도
        fprintf('NASA Blue Marble 이미지 로드 중...\n');
        
        % 여러 소스 시도
        urls = {
            'https://eoimages.gsfc.nasa.gov/images/imagerecords/74000/74393/world.topo.200407.3x5400x2700.jpg',
            'https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73909/world.topo.bathy.200412.3x5400x2700.jpg',
            'https://eoimages.gsfc.nasa.gov/images/imagerecords/57000/57730/land_ocean_ice_2048.jpg'
        };
        
        earth_img = [];
        for i = 1:length(urls)
            try
                earth_img = imread(urls{i});
                fprintf('이미지 로드 성공!\n');
                break;
            catch
                continue;
            end
        end
        
        if isempty(earth_img)
            error('온라인 이미지 로드 실패');
        end
        
        % 이미지 크기 조정 (메모리 효율성을 위해)
        if size(earth_img, 1) > 2000 || size(earth_img, 2) > 4000
            earth_img = imresize(earth_img, [1000 2000]);
        end
        
        % 이미지가 흑백인 경우 RGB로 변환
        if size(earth_img, 3) == 1
            earth_img = repmat(earth_img, [1, 1, 3]);
        end
        
        % 텍스처 매핑으로 지구 그리기
        earth_surf = surf(X, Y, Z, ...
                         'CData', earth_img, ...
                         'FaceColor', 'texturemap', ...
                         'EdgeColor', 'none', ...
                         'FaceLighting', 'gouraud', ...
                         'AmbientStrength', 0.3, ...
                         'DiffuseStrength', 0.8, ...
                         'SpecularStrength', 0.2);
                         
    catch ME
        fprintf('NASA 이미지 로드 실패: %s\n', ME.message);
        fprintf('로컬 파일 사용을 시도합니다...\n');
        
        % 로컬 파일 시도
        local_files = {'earth_texture.jpg', 'blue_marble.jpg', 'earth.jpg'};
        loaded = false;
        
        for i = 1:length(local_files)
            if exist(local_files{i}, 'file')
                try
                    earth_img = imread(local_files{i});
                    fprintf('로컬 이미지 "%s" 로드 성공!\n', local_files{i});
                    
                    % 텍스처 매핑
                    earth_surf = surf(X, Y, Z, ...
                                     'CData', earth_img, ...
                                     'FaceColor', 'texturemap', ...
                                     'EdgeColor', 'none', ...
                                     'FaceLighting', 'gouraud', ...
                                     'AmbientStrength', 0.3, ...
                                     'DiffuseStrength', 0.8, ...
                                     'SpecularStrength', 0.2);
                    loaded = true;
                    break;
                catch
                    continue;
                end
            end
        end
        
        if ~loaded
            fprintf('NASA 이미지를 사용할 수 없습니다. 고품질 대체 텍스처를 생성합니다.\n');
            
            % 고품질 대체 텍스처 생성
            C = generateRealisticEarthTexture(lon, lat);
            
            earth_surf = surf(X, Y, Z, C, ...
                             'EdgeColor', 'none', ...
                             'FaceColor', 'interp', ...
                             'FaceLighting', 'gouraud', ...
                             'AmbientStrength', 0.3, ...
                             'DiffuseStrength', 0.8, ...
                             'SpecularStrength', 0.2);
        end
    end
    
    % 대기층 효과
    atmosphere = surf(X*1.015, Y*1.015, Z*1.015, ...
                     ones(size(X)), ...
                     'EdgeColor', 'none', ...
                     'FaceColor', [0.7 0.85 1], ...
                     'FaceAlpha', 0.06, ...
                     'FaceLighting', 'none');
    
    % 선택적: 주요 위도/경도선
    hold on;
    
    % 적도 (선택적)
    % theta = linspace(0, 2*pi, 100);
    % x_eq = R_earth/1e6 * cos(theta);
    % y_eq = R_earth/1e6 * sin(theta);
    % z_eq = zeros(size(theta));
    % plot3(x_eq, y_eq, z_eq, 'Color', [1 1 1 0.2], 'LineWidth', 0.5);
end

function [r_chief_eci, r_deputy_eci] = computeECITrajectories(X_lvlh, time_vec, chief, params)
    % LVLH 좌표를 ECI 좌표로 변환
    
    N = length(time_vec);
    r_chief_eci = zeros(3, N);
    r_deputy_eci = zeros(3, N);
    
    for k = 1:N
        M = chief.M + chief.n * time_vec(k);
        E = solveKeplerEquation(M, chief.e);
        nu = 2 * atan2(sqrt(1 + chief.e) * sin(E/2), sqrt(1 - chief.e) * cos(E/2));
        
        [r_chief, v_chief] = oe2eci(chief.a, chief.e, chief.i, ...
                                    chief.RAAN, chief.w, nu, params.mu);
        r_chief_eci(:,k) = r_chief;
        
        R_lvlh2eci = getLVLH2ECIMatrix(r_chief, v_chief);
        r_rel_lvlh = X_lvlh([1,3,5], k);
        r_deputy_eci(:,k) = r_chief + R_lvlh2eci * r_rel_lvlh;
    end
end

function E = solveKeplerEquation(M, e)
    % 케플러 방정식 풀이
    E = M;
    tol = 1e-12;
    max_iter = 50;
    
    for i = 1:max_iter
        f = E - e*sin(E) - M;
        fp = 1 - e*cos(E);
        E_new = E - f/fp;
        
        if abs(E_new - E) < tol
            break;
        end
        E = E_new;
    end
end

function [r, v] = oe2eci(a, e, i, RAAN, w, nu, mu)
    % 궤도 요소를 ECI 위치/속도로 변환
    
    p = a * (1 - e^2);
    r_mag = p / (1 + e * cos(nu));
    
    r_orbit = r_mag * [cos(nu); sin(nu); 0];
    
    h = sqrt(mu * p);
    v_orbit = (mu/h) * [-sin(nu); e + cos(nu); 0];
    
    R3_RAAN = [cos(RAAN), -sin(RAAN), 0; sin(RAAN), cos(RAAN), 0; 0, 0, 1];
    R1_i = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)];
    R3_w = [cos(w), -sin(w), 0; sin(w), cos(w), 0; 0, 0, 1];
    
    R_orbit2eci = R3_RAAN * R1_i * R3_w;
    
    r = R_orbit2eci * r_orbit;
    v = R_orbit2eci * v_orbit;
end

function R = getLVLH2ECIMatrix(r, v)
    % LVLH to ECI 변환 행렬
    
    z_hat = -r / norm(r);
    h = cross(r, v);
    y_hat = -h / norm(h);
    x_hat = cross(y_hat, z_hat);
    
    R = [x_hat'; y_hat'; z_hat']';
end

function addStars(ax)
    % 배경 별 추가
    n_stars = 500;
    
    theta = 2*pi*rand(n_stars, 1);
    phi = acos(2*rand(n_stars, 1) - 1);
    r = 15000;
    
    x_stars = r * sin(phi) .* cos(theta);
    y_stars = r * sin(phi) .* sin(theta);
    z_stars = r * cos(phi);
    
    brightness = 0.3 + 0.7*rand(n_stars, 1);
    
    scatter3(ax, x_stars/1e3, y_stars/1e3, z_stars/1e3, ...
             brightness*3, [1 1 1], 'filled', 'MarkerFaceAlpha', 0.8);
end