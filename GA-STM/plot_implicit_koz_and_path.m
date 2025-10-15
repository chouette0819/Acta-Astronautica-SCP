%% Plot all implicit KOZ surfaces F(x,y,z)-1=0 and overlay trajectory
clear; clc; close all;

% Load SCP and (optionally) QP results
thisdir = fileparts(mfilename('fullpath'));
scp_mat = fullfile(thisdir,'scp_last.mat');
qp_mat  = fullfile(thisdir,'qp_last.mat');
has_scp = exist(scp_mat,'file')~=0; has_qp = exist(qp_mat,'file')~=0;
if ~has_scp
    error('scp_last.mat not found. Run SCP_KOZ_QP.m first.');
end
S = load(scp_mat); X_scp = S.X; t = S.time_vec; %#ok<NASGU>
if has_qp
    Q = load(qp_mat); X_qp = Q.X_opt; %#ok<NASGU>
end

% Determine plotting bounds: fixed cube [-100,100] m on each axis
P = X_scp([1,3,5],:); % [x;y;z]
lo = [-100; -100; -100]; % [m]
hi = [ 100;  100;  100]; % [m]

% Base grid settings (higher resolution)
res = 100; % increase for finer surfaces
nx=res; ny=res; nz=res;

% Discover implicit surface expressions
expr_dir = 'C:\Users\98kim\Desktop\Acta-Astronautica\Funcs_ISS_expr';
if ~exist(expr_dir,'dir')
    error('Expr dir not found: %s', expr_dir);
end
files = dir(fullfile(expr_dir,'*_*_expr.txt'));
max_surfaces = numel(files);
if isempty(files)
    error('No *_expr.txt found in %s', expr_dir);
end

% Figure
fig = figure('Position',[80 80 1200 900]); hold on; grid on;
set(gcf,'Renderer','opengl');
colormap(parula);
colors = lines(max_surfaces);

% Plot all implicit surfaces
for i=1:max_surfaces
    expr_path = fullfile(expr_dir, files(i).name);
    Fh = make_F_from_expr(expr_path);
    % Adaptive grid per-surface to ensure F-1=0 is present if nearby
    lo_i = lo; hi_i = hi; drawn = false;
    for expand_it = 1:3 % try up to 3 expansions
        xv = linspace(lo_i(1),hi_i(1),nx);
        yv = linspace(lo_i(2),hi_i(2),ny);
        zv = linspace(lo_i(3),hi_i(3),nz);
        [Xg,Yg,Zg] = meshgrid(xv,yv,zv);
        V = arrayfun(@(x,y,z) Fh(x,y,z), Xg, Yg, Zg) - 1; % F-1 (m)
        vmin = min(V(:)); vmax = max(V(:));
        fprintf('Surface %d (%s): min=%.2e max=%.2e\n', i, files(i).name, vmin, vmax);
        if vmin<=0 && vmax>=0
            try
                p = patch(isosurface(Xg/1e3,Yg/1e3,Zg/1e3,V,0));
                c = colors(1+mod(i-1,size(colors,1)),:);
                set(p,'FaceColor',c,'EdgeColor','none','FaceAlpha',0.35);
                isonormals(Xg/1e3,Yg/1e3,Zg/1e3,V,p);
            drawn = true; break;
            catch
                % fall through to expand
            end
        end
        % expand box by 2x around center and retry
        center = 0.5*(lo_i+hi_i);
        half = 0.5*(hi_i-lo_i);
        lo_i = center - 2*half; hi_i = center + 2*half;
    end
    if ~drawn
        fprintf('  (no F=1 isosurface found near path bbox for %s)\n', files(i).name);
    end
end

% Overlay trajectories
plot3(P(1,:)/1e3, P(2,:)/1e3, P(3,:)/1e3,'r-','LineWidth',2);
if has_qp
    plot3(X_qp(1,:)/1e3, X_qp(3,:)/1e3, X_qp(5,:)/1e3,'b-','LineWidth',2);
end

% Formatting
axis equal; view(45,30);
xlim([-0.1 0.1]); ylim([-0.1 0.1]); zlim([-0.1 0.1]); % +/-100 m in km
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
legend_entries = {'Implicit surfaces','SCP path'};
if has_qp, legend_entries{end+1} = 'QP path'; end %#ok<AGROW>
legend(legend_entries,'Location','best');
title('All implicit KOZ surfaces F-1=0 with trajectory');
lighting gouraud; camlight headlight; material dull;

% Save
print(fig, fullfile(thisdir,'koz_surfaces_with_path.png'), '-dpng','-r150');

%% Helpers
function Fh = make_F_from_expr(expr_path)
    txt = strtrim(fileread(expr_path));
    txt = regexprep(txt, '\*\*', '.^'); % Python power to MATLAB
    fstr = sprintf('@(x,y,z) (%s)', txt);
    Fh = str2func(fstr);
end
