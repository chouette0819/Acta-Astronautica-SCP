%% Compare plain QP vs SCP (KOZ) trajectories on one figure
clear; clc; close all;
thisdir = fileparts(mfilename('fullpath'));

% Load precomputed results saved by the scripts
scp_mat = fullfile(thisdir,'scp_last.mat');
qp_mat  = fullfile(thisdir,'qp_last.mat');
if ~exist(scp_mat,'file') || ~exist(qp_mat,'file')
    error('Missing scp_last.mat or qp_last.mat. Run SCP_KOZ_QP.m and OptimizationBaselineQP.m first.');
end
S = load(scp_mat); Q = load(qp_mat);
X_scp = S.X; U_scp = S.U; koz_scp = S.koz;
X_qp  = Q.X_opt;

% Plot
fig = figure('Position',[80 80 900 700]);
hold on; grid on;
plot3(X_qp(1,:)/1e3, X_qp(3,:)/1e3, X_qp(5,:)/1e3,'b-','LineWidth',2);
plot3(X_scp(1,:)/1e3, X_scp(3,:)/1e3, X_scp(5,:)/1e3,'r-','LineWidth',2);
% Do not draw ellipsoidal KOZ here; only trajectories
axis equal; view(45,30);
xlim([-0.1 0.1]); ylim([-0.1 0.1]); zlim([-0.1 0.1]); % +/-100 m in km
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend({'QP','SCP (KOZ)'},'Location','best');
title('Trajectory Comparison: QP vs SCP (KOZ)');

% Save
print(fig,'compare_qp_scp.png','-dpng','-r150');
