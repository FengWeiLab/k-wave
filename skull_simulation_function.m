function result = skull_simulation_function(cfg, opts)
%SKULL_SIMULATION_FUNCTION 颅骨超声相控阵 k-Wave 仿真主函数
%
%  文件 2/3：skull_simulation_function.m
%
%  本文件负责 simulate_32el_15MHz_curved_skull.m 的全部仿真功能：
%    1) 按配置建立 OACA 处理前 / 处理后的三层曲面介质模型（含缺陷）；
%    2) 导出模型声速、密度分布 MATLAB Fig 图（2_skull_model，按模板版式）；
%    3) 逐阵元发射、全阵接收的 k-Wave 二维仿真（FMC 全矩阵采集）；
%    4) 保存 OACA 处理前后的相控阵数据 mat 文件（3_signal，变量与
%       simdata_32el_15MHz_2um_curved_skull_sphere_*.mat 完全一致）；
%    5) 录制并保存 OACA 处理前后的声场传播动画 avi（5_video）。
%
%  用法：
%    result = skull_simulation_function(CFG)                  % 按 CFG.precheck 自动选择
%    result = skull_simulation_function(CFG, struct('mode', 'model'))     % 只建模出图
%    result = skull_simulation_function(CFG, struct('mode', 'simulate'))  % 完整仿真
%
%  输入：
%    cfg  ：skull_simulation_model.m 中定义的配置结构体
%    opts ：可选控制项
%             .mode      'auto'（默认）| 'model' | 'simulate'
%             .cases     需要运行的 case tag 元胞数组（默认 cfg.run_cases）
%             .overwrite 是否覆盖已有仿真数据（默认 false，存在则跳过）
%
%  输出：
%    result.folder   4 个输出目录的绝对路径
%    result.cases    每个处理状态的产物清单
%    result.setup    网格/阵列等公共设置（便于其它脚本复用）
%
%  说明：
%    * 每个阵元用 1 个网格点作为点源与点接收，32 阵元共 32 次发射；
%    * 时间步 dt = dx / min(c) / dt_div，采样率 fs = 1/dt；
%    * 长仿真（单套配置约数十小时）支持逐阵元断点续算。
%
%  参考：simulate_32el_15MHz_curved_skull.m
% =========================================================================

%% ---------------- 0. 输入解析 ----------------
if nargin < 1 || isempty(cfg)
    error('skull_simulation_function:no_cfg', ...
        '必须提供配置结构体 CFG，请先运行 skull_simulation_model.m。');
end
if nargin < 2 || isempty(opts)
    opts = struct();
end

cfg = local_normalize(cfg);

if ~isfield(opts, 'mode') || isempty(opts.mode)
    if cfg.precheck == 1
        opts.mode = 'simulate';
    else
        opts.mode = 'model';
    end
end
if ~isfield(opts, 'cases') || isempty(opts.cases)
    opts.cases = cfg.run_cases;
end
if ~isfield(opts, 'overwrite') || isempty(opts.overwrite)
    opts.overwrite = false;
end
if ~any(strcmp(opts.mode, {'model', 'simulate'}))
    error('skull_simulation_function:bad_mode', ...
        'opts.mode 只能是 ''model'' 或 ''simulate''，当前为 ''%s''。', opts.mode);
end

do_simulate = strcmp(opts.mode, 'simulate');

%% ---------------- 1. k-Wave 工具箱路径 ----------------
kpath = cfg.kwave_dir;
if ~exist(kpath, 'dir')
    error('skull_simulation_function:no_kwave', '未找到 k-Wave 工具箱：%s', kpath);
end
addpath(kpath);

%% ---------------- 2. 输出目录 ----------------
F = cfg.folder;
folder = struct();
folder.model = fullfile(cfg.root_dir, F.model);
folder.signal = fullfile(cfg.root_dir, F.signal);
folder.image = fullfile(cfg.root_dir, F.image);
folder.video = fullfile(cfg.root_dir, F.video);
fn = fieldnames(folder);
for i = 1:numel(fn)
    if ~exist(folder.(fn{i}), 'dir')
        mkdir(folder.(fn{i}));
    end
end
result = struct();
result.folder = folder;

%% ---------------- 3. 公共设置（网格、阵列、计算精度）----------------
S = build_setup(cfg);
S.data_cast = detect_data_cast(cfg.sim.use_gpu);

fprintf('\n---------------- 仿真设置 ----------------\n');
fprintf('计算域：Nx = %d（深度 %.3f mm），Ny = %d（横向 %.3f mm）\n', ...
    S.Nx, S.Nx * cfg.grid.dx * 1e3, S.Ny, S.Ny * cfg.grid.dx * 1e3);
fprintf('阵列  ：%d 阵元，间距 %.1f um，孔径 %.2f mm，深度 %.1f um\n', ...
    cfg.array.n_elements, cfg.array.pitch * 1e6, ...
    (cfg.array.n_elements - 1) * cfg.array.pitch * 1e3, S.array_depth * 1e6);
fprintf('精度  ：%s\n', S.data_cast);
fprintf('模式  ：%s\n', opts.mode);
fprintf('------------------------------------------\n');
result.setup = S;

%% ---------------- 4. 逐状态运行 ----------------
case_summary = struct('tag', {}, 'label', {}, 'full_data', {}, ...
    'movie', {}, 'model_figures', {}, 'elapsed', {});

t_all = tic;
for ic = 1:numel(opts.cases)
    tag = opts.cases{ic};
    k = find(strcmp({cfg.cases.tag}, tag), 1);
    if isempty(k)
        warning('skull_simulation_function:unknown_case', ...
            '未定义的配置标签 %s，已跳过。', tag);
        continue;
    end
    case_def = cfg.cases(k);

    fprintf('\n================ [%d/%d] %s ================\n', ...
        ic, numel(opts.cases), case_def.label);
    t_case = tic;

    % ---- 4.1 建立介质模型（含缺陷）----
    [medium, geom] = build_medium(cfg, case_def, S);
    print_medium_info(cfg, case_def, medium, geom, S);

    % ---- 4.2 导出模型分布图（2_skull_model）----
    fig_files = {};
    if cfg.output.model_figures
        fig_files = plot_model_maps(cfg, case_def, medium, geom, S, folder.model);
    end

    full_data_file = '';
    movie_file = '';

    if do_simulate
        % ---- 4.3 逐阵元发射 / 全阵接收（FMC）----
        % 相控阵数据与 simdata_*.mat 同格式，存在 3_signal
        full_data_file = fullfile(folder.signal, ...
            sprintf('%s_full_data.mat', case_def.short));
        if exist(full_data_file, 'file') && ~opts.overwrite
            fprintf('  已存在相控阵数据，跳过仿真：%s\n', full_data_file);
            fprintf('  （需要重算请设置 opts.overwrite = true，或删除该文件）\n');
        else
            [~, sim_config, input_signal, t_axis] = run_fmc( ...
                cfg, case_def, S, medium, geom, full_data_file);
            % ---- 4.4 声场传播动画 ----
            if cfg.save_animation
                movie_file = run_animation(cfg, case_def, S, medium, ...
                    sim_config, input_signal, t_axis, folder.video);
            end
        end
    end

    elapsed = toc(t_case);
    fprintf('  %s 完成，用时 %.1f s\n', case_def.label, elapsed);

    case_summary(end + 1) = struct('tag', case_def.tag, ...
        'label', case_def.label, 'full_data', full_data_file, ...
        'movie', movie_file, 'model_figures', {fig_files}, ...
        'elapsed', elapsed); %#ok<SAGROW>
end

result.cases = case_summary;
result.elapsed = toc(t_all);
result.mode = opts.mode;

%% ---------------- 5. 收尾汇总 ----------------
fprintf('\n================ 仿真阶段汇总（%.1f s）================\n', result.elapsed);
for i = 1:numel(case_summary)
    c = case_summary(i);
    fprintf('[%s]\n', c.label);
    for j = 1:numel(c.model_figures)
        fprintf('    %s\n', c.model_figures{j});
    end
    if ~isempty(c.full_data)
        fprintf('    %s\n', c.full_data);
    end
    if ~isempty(c.movie)
        fprintf('    %s\n', c.movie);
    end
end
fprintf('======================================================\n');

end % ==================== 主函数结束 ====================


%% =====================================================================
%% 配置规范化
%% =====================================================================
function cfg = local_normalize(cfg)
%LOCAL_NORMALIZE 补齐缺省字段，保证老配置也能运行

if ~isfield(cfg, 'precheck') || isempty(cfg.precheck)
    cfg.precheck = 0;
end
if ~isfield(cfg, 'save_animation') || isempty(cfg.save_animation)
    cfg.save_animation = 0;
end
if ~isfield(cfg, 'run_cases') || isempty(cfg.run_cases)
    cfg.run_cases = {cfg.cases.tag};
end
if ~isfield(cfg, 'folder') || isempty(cfg.folder)
    cfg.folder = struct('model', '2_skull_model', 'signal', '3_signal', ...
        'image', '4_Image', 'video', '5_video');
end
if ~isfield(cfg, 'code_dir') || isempty(cfg.code_dir)
    cfg.code_dir = fileparts(mfilename('fullpath'));
end
if ~isfield(cfg, 'root_dir') || isempty(cfg.root_dir)
    cfg.root_dir = fileparts(cfg.code_dir);         % 代码目录的上一级
    if isempty(cfg.root_dir) || ~exist(cfg.root_dir, 'dir')
        cfg.root_dir = cfg.code_dir;
    end
end
if ~isfield(cfg, 'sim') || isempty(cfg.sim)
    cfg.sim = struct('use_gpu', true, 'save_checkpoint', true, ...
        'checkpoint_every', 8, 'resume', true, 'save_full_data', true, ...
        'verbose', true);
else
    d = struct('use_gpu', true, 'save_checkpoint', true, ...
        'checkpoint_every', 8, 'resume', true, 'save_full_data', true, ...
        'verbose', true);
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(cfg.sim, f{i}) || isempty(cfg.sim.(f{i}))
            cfg.sim.(f{i}) = d.(f{i});
        end
    end
end
if ~isfield(cfg, 'output') || isempty(cfg.output)
    cfg.output = struct('model_figures', true, 'export_png', true, ...
        'close_figures', true);
else
    d = struct('model_figures', true, 'export_png', true, ...
        'close_figures', true);
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(cfg.output, f{i}) || isempty(cfg.output.(f{i}))
            cfg.output.(f{i}) = d.(f{i});
        end
    end
end
% 模型分布图版式（与 2_skull_model\OACA_Speed.fig 模板一致）
d_mp = struct('clim_speed', [1400 3400], 'clim_density', [900 2000], ...
    'clim_impedance', [1.4e6 6.2e6], 'font_name', 'Times New Roman', ...
    'font_size', 12, 'figure_size', [277.6 404], 'draw_impedance', false);
if ~isfield(cfg, 'model_plot') || isempty(cfg.model_plot)
    cfg.model_plot = d_mp;
else
    f = fieldnames(d_mp);
    for i = 1:numel(f)
        if ~isfield(cfg.model_plot, f{i}) || isempty(cfg.model_plot.(f{i}))
            cfg.model_plot.(f{i}) = d_mp.(f{i});
        end
    end
end
if ~isfield(cfg, 'animation') || isempty(cfg.animation)
    cfg.animation = struct('aperture', 'center', 'plot_freq', 60, ...
        'frame_rate', 25, 'plot_scale', 'auto', 'profile', 'Motion JPEG AVI', ...
        'figure_position', [100 50 560 860]);
else
    d = struct('aperture', 'center', 'plot_freq', 60, 'frame_rate', 25, ...
        'plot_scale', 'auto', 'profile', 'Motion JPEG AVI', ...
        'figure_position', [100 50 560 860]);
    f = fieldnames(d);
    for i = 1:numel(f)
        if ~isfield(cfg.animation, f{i}) || isempty(cfg.animation.(f{i}))
            cfg.animation.(f{i}) = d.(f{i});
        end
    end
end
if ~isfield(cfg, 'defects') || isempty(cfg.defects)
    cfg.defects = struct('depth', 5.0e-3, 'lateral', 0, 'diameter', 0.1e-3, ...
        'c', 2000, 'rho', 1500, 'label', '散射体');
end
if ~isfield(cfg, 'kwave_dir') || isempty(cfg.kwave_dir)
    cfg.kwave_dir = fullfile(cfg.code_dir, 'k-wave-toolbox-version-1.4', 'k-Wave');
end
end


%% =====================================================================
%% 网格与阵列
%% =====================================================================
function S = build_setup(cfg)
%BUILD_SETUP 建立计算网格坐标、阵列位置与采样参数

P = cfg.grid;
Nx = round((cfg.geometry.gel_thickness + cfg.geometry.skull_thickness + ...
    cfg.geometry.brain_thickness) / P.dx);
Ny = round(P.lateral_width / P.dx);

S.Nx = Nx;
S.Ny = Ny;
S.dx = P.dx;
S.dy = P.dx;
S.x_axis = (0:Nx - 1) * P.dx;                       % 深度坐标 [m]
S.y_axis = ((1:Ny) - (Ny + 1) / 2) * P.dx;          % 横向坐标 [m]

% 阵列：每个阵元一个网格点，放在顶面下方一个网格点处，避免源点落在边界上
S.array_depth_idx = 2;
array_y = ((1:cfg.array.n_elements) - (cfg.array.n_elements + 1) / 2) * cfg.array.pitch;
array_y_idx = round((array_y - S.y_axis(1)) / S.dy) + 1;
S.array_y_requested = array_y;
S.array_y_idx = max(1, min(Ny, array_y_idx));
S.array_y_actual = S.y_axis(S.array_y_idx);
S.array_depth = S.x_axis(S.array_depth_idx);

if numel(unique(S.array_y_idx)) ~= cfg.array.n_elements
    error('skull_simulation_function:array_overlap', ...
        ['阵元在网格上发生重叠：%d 个阵元仅占 %d 个网格点，' ...
        '请减小 pitch 或阵元数。'], cfg.array.n_elements, ...
        numel(unique(S.array_y_idx)));
end
end


%% =====================================================================
%% 计算精度 / GPU
%% =====================================================================
function data_cast = detect_data_cast(use_gpu)
%DETECT_DATA_CAST 有可用 GPU 时返回 'gpuArray-single'，否则返回 'single'

data_cast = 'single';
if ~use_gpu
    fprintf('按配置使用 CPU 计算\n');
    return
end
try
    has_func = (exist('gpuDeviceCount', 'builtin') > 0) || ...
        (exist('gpuDeviceCount', 'file') > 0);
    if has_func && gpuDeviceCount > 0
        data_cast = 'gpuArray-single';
        g = gpuDevice;
        fprintf('检测到 GPU：%s（%.1f GB）\n', g.Name, g.TotalMemory / 2^30);
    else
        fprintf('未检测到可用 GPU，使用 CPU 计算\n');
    end
catch
    data_cast = 'single';
    fprintf('GPU 检测失败，使用 CPU 计算\n');
end
end


%% =====================================================================
%% 介质建模
%% =====================================================================
function [medium, geom] = build_medium(cfg, case_def, S)
%BUILD_MEDIUM 建立三层曲面介质，并在脑层内嵌入缺陷（散射体）
%
%   球面矢高 sag(R, y) = R - sqrt(R^2 - y^2)：
%     凝胶下表面 = 凝胶厚度 + sag(R1)
%     颅骨下表面 = 凝胶厚度 + 颅骨厚度 + sag(R2)
%   默认填充脑层；计算域底边即脑层下表面（平面）。

P = cfg.geometry;
x_axis = S.x_axis;
y_axis = S.y_axis;
Nx = numel(x_axis);
Ny = numel(y_axis);

sag_gel_skull = P.radius_gel_skull - ...
    sqrt(max(P.radius_gel_skull^2 - y_axis.^2, 0));
sag_skull_brain = P.radius_skull_brain - ...
    sqrt(max(P.radius_skull_brain^2 - y_axis.^2, 0));

geom.z_gel_bottom = P.gel_thickness + sag_gel_skull;
geom.z_skull_bottom = P.gel_thickness + P.skull_thickness + sag_skull_brain;
geom.z_brain_bottom = P.gel_thickness + P.skull_thickness + P.brain_thickness;

% 默认填充脑层
medium.sound_speed = single(cfg.media.c_brain * ones(Nx, Ny, 'single'));
medium.density = single(cfg.media.rho_brain * ones(Nx, Ny, 'single'));

for iy = 1:Ny
    gel_end = max(1, min(Nx, floor(geom.z_gel_bottom(iy) / S.dx) + 1));
    skull_end = max(gel_end + 1, min(Nx, floor(geom.z_skull_bottom(iy) / S.dx) + 1));

    medium.sound_speed(1:gel_end, iy) = single(cfg.media.c_gel);
    medium.density(1:gel_end, iy) = single(cfg.media.rho_gel);

    if skull_end > gel_end
        medium.sound_speed(gel_end + 1:skull_end, iy) = single(case_def.c_skull);
        medium.density(gel_end + 1:skull_end, iy) = single(case_def.rho_skull);
    end
end

% ---- 缺陷（球形散射体，二维为圆截面）----
n_def = numel(cfg.defects);
geom.defects = struct('label', {}, 'depth', {}, 'lateral', {}, 'diameter', {}, ...
    'c', {}, 'rho', {}, 'ix', {}, 'iy', {}, 'radius_px', {}, 'npix', {});

for id = 1:n_def
    d = cfg.defects(id);
    r = 0.5 * d.diameter;
    r_px = r / S.dx;

    % 缺陷中心在连续网格坐标下的位置（行=深度，列=横向）
    ix = d.depth / S.dx + 1;
    iy = (d.lateral - y_axis(1)) / S.dy + 1;

    pad = ceil(r_px) + 2;
    rows = max(1, floor(ix - pad)):min(Nx, ceil(ix + pad));
    cols = max(1, floor(iy - pad)):min(Ny, ceil(iy + pad));
    [Rsub, Csub] = ndgrid(rows, cols);
    disk = (Rsub - ix).^2 + (Csub - iy).^2 <= r_px^2;

    c_sub = medium.sound_speed(rows, cols);
    rho_sub = medium.density(rows, cols);
    c_sub(disk) = single(d.c);
    rho_sub(disk) = single(d.rho);
    medium.sound_speed(rows, cols) = c_sub;
    medium.density(rows, cols) = rho_sub;

    geom.defects(id) = struct('label', d.label, 'depth', d.depth, ...
        'lateral', d.lateral, 'diameter', d.diameter, 'c', d.c, 'rho', d.rho, ...
        'ix', ix, 'iy', iy, 'radius_px', r_px, 'npix', nnz(disk));

    % ---- 位置检查：缺陷应完整落在脑层内部 ----
    iy0 = max(1, min(Ny, round(iy)));
    if d.depth - r <= geom.z_skull_bottom(iy0)
        warning('skull_simulation_function:defect_overlap', ...
            '缺陷 [%d] %s 与颅骨层重叠：上缘 %.3f mm，颅骨下表面 %.3f mm', ...
            id, d.label, (d.depth - r) * 1e3, geom.z_skull_bottom(iy0) * 1e3);
    end
    if d.depth + r >= geom.z_brain_bottom
        warning('skull_simulation_function:defect_outside', ...
            '缺陷 [%d] %s 超出脑层下表面：下缘 %.3f mm，脑层下表面 %.3f mm', ...
            id, d.label, (d.depth + r) * 1e3, geom.z_brain_bottom * 1e3);
    end
    if d.lateral - r < y_axis(1) || d.lateral + r > y_axis(end)
        warning('skull_simulation_function:defect_outside', ...
            '缺陷 [%d] %s 超出横向计算范围', id, d.label);
    end
    if geom.defects(id).npix < 1
        warning('skull_simulation_function:defect_empty', ...
            '缺陷 [%d] %s 未落到任何网格点上，请检查直径/位置', id, d.label);
    end
end
end


function print_medium_info(cfg, case_def, medium, geom, S)
%PRINT_MEDIUM_INFO 打印介质与缺陷信息（含合理性检查）

fprintf('  介质：凝胶 %.0f/%.0f，颅骨 %.1f/%.0f，脑 %.0f/%.0f  (c/rho)\n', ...
    cfg.media.c_gel, cfg.media.rho_gel, case_def.c_skull, ...
    case_def.rho_skull, cfg.media.c_brain, cfg.media.rho_brain);
fprintf('  层界面（中心处）：凝胶/颅骨 %.3f mm，颅骨/脑 %.3f mm，脑底 %.3f mm\n', ...
    geom.z_gel_bottom(round((S.Ny + 1) / 2)) * 1e3, ...
    geom.z_skull_bottom(round((S.Ny + 1) / 2)) * 1e3, ...
    geom.z_brain_bottom * 1e3);
for id = 1:numel(geom.defects)
    dd = geom.defects(id);
    fprintf('  缺陷 [%d] %s：深度 %.3f mm，横向 %.3f mm，直径 %.3f mm，占 %d 个网格点\n', ...
        id, dd.label, dd.depth * 1e3, dd.lateral * 1e3, ...
        dd.diameter * 1e3, dd.npix);
end
fprintf('  声速范围 %.0f ~ %.0f m/s，密度范围 %.0f ~ %.0f kg/m^3\n', ...
    min(medium.sound_speed(:)), max(medium.sound_speed(:)), ...
    min(medium.density(:)), max(medium.density(:)));
fprintf('  网格：Nx = %d，Ny = %d\n', S.Nx, S.Ny);
end


%% =====================================================================
%% 模型分布图（文件夹 2_skull_model）
%% =====================================================================
function files = plot_model_maps(cfg, case_def, medium, geom, S, out_dir)
%PLOT_MODEL_MAPS 导出声速 / 密度分布 MATLAB Fig 图
%
%   版式按照模板 2_skull_model\OACA_Speed.fig：
%     图窗 277.6 x 404 px（CFG.model_plot.figure_size），
%     Times New Roman 12，Box on，TickDir in，
%     标题 "Speed before/after OACA"、"Density before/after OACA"，
%     坐标 "x / mm"、"z / mm"，
%     色标单位用图窗文本框放在色标正下方（"m/s"、"kg/m^3"）。
%   坐标区位置交给 MATLAB 自动布局（模板里存的 [0 ...] 会把纵轴刻度裁到
%   图窗外面），其余字号、色标、标注位置均按模板。
%   处理前 / 处理后的色标范围固定一致（见 CFG.model_plot）：
%     声速 1400~3400 m/s，密度 900~2000 kg/m^3。
%   声阻抗按最新要求不再出图（CFG.model_plot.draw_impedance = false）。

mp = cfg.model_plot;
c = double(medium.sound_speed);
rho = double(medium.density);

% 处理前 / 处理后（用于英文标题 "Speed before OACA" / "Speed after OACA"）
if isfield(case_def, 'en') && ~isempty(case_def.en)
    when = case_def.en;
elseif strcmpi(case_def.short, 'original')
    when = 'before';
else
    when = 'after';
end

specs = { ...
    'Speed',   c,   'Speed',   'm/s',     mp.clim_speed; ...
    'Density', rho, 'Density', 'kg/m^3',  mp.clim_density};

if isfield(mp, 'draw_impedance') && mp.draw_impedance
    specs = [specs; {'Impedance', c .* rho, 'Impedance', ...
        'kg/(m^2 s)',  mp.clim_impedance}]; %#ok<AGROW>
end

files = cell(1, size(specs, 1));
for i = 1:size(specs, 1)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', ...
        [100 100 mp.figure_size]);
    ax = axes(fig);

    ttl = sprintf('%s %s OACA', specs{i, 3}, when);
    draw_map(ax, cfg, S, geom, specs{i, 2}, ttl, specs{i, 5}, specs{i, 4});

    base = fullfile(out_dir, sprintf('%s_%s', case_def.short, specs{i, 1}));
    saved = save_figure(fig, base, cfg.output.export_png, cfg.output.close_figures, 300);
    files{i} = saved{1};
    fprintf('  已保存模型图：%s.fig\n', base);
end
end


function draw_map(ax, cfg, S, geom, data, ttl, clim_range, unit)
%DRAW_MAP 绘制一张材料分布图（版式按 2_skull_model\OACA_Speed.fig 模板）
%
%   显示方向：阵列在上方（深度 0），大脑在下方，深度向下增大。
%   注意：这里不使用 colorbar 的 Label 属性（在图窗更新遍历中访问
%   cb.Label 会触发 "尝试在更新遍历期间修改树"），单位改用图窗文本框，
%   放在色标下方，与模板一致。

mp = cfg.model_plot;

imagesc(ax, S.y_axis * 1e3, S.x_axis * 1e3, data);
set(ax, 'YDir', 'reverse');
% 顶部留一点空白，避免阵元标记贴边被裁掉（与模板 YLim 一致）
xlim(ax, [S.y_axis(1) S.y_axis(end)] * 1e3);
ylim(ax, [S.x_axis(1) * 1e3 - 0.06, S.x_axis(end) * 1e3]);
hold(ax, 'on');

% 层界面（凝胶/颅骨 白线，颅骨/大脑 黑线，脑底 品红线）
plot(ax, S.y_axis * 1e3, geom.z_gel_bottom * 1e3, 'w-', 'LineWidth', 1.0);
plot(ax, S.y_axis * 1e3, geom.z_skull_bottom * 1e3, 'k-', 'LineWidth', 1.0);
plot(ax, [S.y_axis(1) S.y_axis(end)] * 1e3, geom.z_brain_bottom * 1e3 * [1 1], ...
    'm-', 'LineWidth', 1.0);

% 32 个阵元（红色点）
plot(ax, S.array_y_actual * 1e3, S.array_depth * 1e3 * ones(size(S.array_y_actual)), ...
    'r.', 'MarkerSize', 12);

% 缺陷（红色圆圈 + 白色说明，与模板同一位置/字号）
th = linspace(0, 2 * pi, 200);
for id = 1:numel(geom.defects)
    d = geom.defects(id);
    r_mm = 0.5 * d.diameter * 1e3;
    plot(ax, d.lateral * 1e3 + r_mm * cos(th), ...
        d.depth * 1e3 + r_mm * sin(th), 'r-', 'LineWidth', 1.2);
    text(ax, d.lateral * 1e3 - 0.529, d.depth * 1e3 + 0.317, ...
        {'   Scatter', sprintf('d = %.2f mm', d.diameter * 1e3)}, ...
        'Color', 'w', 'FontSize', 8, 'FontName', mp.font_name, ...
        'VerticalAlignment', 'middle', 'Interpreter', 'tex');
end

% 坐标区版式（字体/刻度/边框取自模板；位置交给 MATLAB 自动布局，
% 这样纵轴刻度与 "z / mm" 才不会被图窗左边界裁掉）
set(ax, 'FontName', mp.font_name, 'FontSize', mp.font_size, ...
    'TickDir', 'in', 'Box', 'on', 'Layer', 'top');
title(ax, ttl);
xlabel(ax, 'x / mm');
ylabel(ax, 'z / mm');
if ~isempty(clim_range)
    set(ax, 'CLim', clim_range);
end

cb = colorbar(ax);
set(cb, 'FontName', mp.font_name, 'FontSize', mp.font_size - 2);

% 色标单位：图窗文本框，放在色标正下方（与模板的光标位置一致）
cpos = get(cb, 'Position');
fig = ancestor(ax, 'figure');
annotation(fig, 'textbox', [cpos(1), cpos(2) - 0.078, 0.147, 0.064], ...
    'String', unit, 'FontName', mp.font_name, 'FontSize', mp.font_size - 2, ...
    'Color', [0 0 0], 'EdgeColor', 'none', 'BackgroundColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'Interpreter', 'tex', 'FitBoxToText', 'off');
end


function paths = save_figure(fig, base, export_png, close_after, png_res)
%SAVE_FIGURE 保存 .fig，并可选导出同名 .png 预览
%
%   注意：绘图阶段图窗是 hidden 的（批量出图不弹窗），但保存前必须把
%   Visible 属性改成 on，否则用户双击打开的 .fig 窗口是隐藏的、看不见。

if nargin < 4 || isempty(close_after)
    close_after = true;
end
if nargin < 5 || isempty(png_res)
    png_res = 200;
end
was_visible = strcmpi(get(fig, 'Visible'), 'on');
set(fig, 'Visible', 'on');              % 让 .fig 里保存的属性为可见
savefig(fig, [base '.fig']);
if ~was_visible
    set(fig, 'Visible', 'off');         % 保存后立刻隐藏，避免窗口闪烁
end
paths = {[base '.fig']};
if export_png
    try
        exportgraphics(fig, [base '.png'], 'Resolution', png_res);
        paths{end + 1} = [base '.png']; %#ok<AGROW>
    catch ME
        warning('skull_simulation_function:png', ...
            '导出 png 失败（%s），已跳过：%s.png', ME.message, base);
    end
end
if close_after
    close(fig);
end
end


%% =====================================================================
%% k-Wave 选项
%% =====================================================================
function opts = kwave_opts(cfg, data_cast, extra)
%KWAVE_OPTS 组装 kspaceFirstOrder2D 的输入选项

opts = {'PMLInside', false, ...
    'PMLSize', cfg.grid.pml_size, ...
    'PMLAlpha', cfg.grid.pml_alpha, ...
    'PlotLayout', false, ...
    'PlotPML', false, ...
    'PlotSim', false, ...
    'RecordMovie', false, ...
    'DataCast', data_cast};

if nargin >= 3 && ~isempty(extra)
    opts = [opts, extra]; %#ok<AGROW>
end
end


%% =====================================================================
%% FMC 仿真（逐阵元发射 / 全阵接收）
%% =====================================================================
function [rf, sim_config, input_signal, t_axis] = run_fmc(cfg, case_def, S, medium, geom, mat_file)
%RUN_FMC 逐阵元发射、全阵元接收，返回 rf(t, rx, tx)

Nx = S.Nx;
Ny = S.Ny;
n_tx = cfg.array.n_elements;
data_cast = S.data_cast;

t_start = tic;

% ---- 仿真时间步进（不变）：dt = dx / min(c) / dt_div ----
c_min = min([cfg.media.c_gel, cfg.media.c_brain, case_def.c_skull, [cfg.defects.c]]);
dt = cfg.grid.dx / c_min / cfg.grid.dt_div;

kgrid = kWaveGrid(Nx, cfg.grid.dx, Ny, S.dy);
kgrid.dt = dt;
kgrid.Nt = ceil(cfg.grid.t_end / dt) + 1;
t_native = (0:kgrid.Nt - 1) * dt;
fs_native = 1 / dt;

% ---- 导出采样：仿真步进不变，保存前把数据重采样到 fs_out ----
fs_out = fs_native;
if isfield(cfg.grid, 'fs_out') && ~isempty(cfg.grid.fs_out)
    fs_out = min(cfg.grid.fs_out, fs_native);
end
n_out = floor(t_native(end) * fs_out) + 1;
t_axis = (0:n_out - 1).' / fs_out;

% 网格坐标（随相控阵数据一起保存，便于后处理直接取用）
x_axis = S.x_axis; %#ok<NASGU>
y_axis = S.y_axis; %#ok<NASGU>

% ---- 激励信号（在仿真步进上生成，不受导出采样率影响）----
input_signal = single(cfg.array.source_pressure * toneBurst(fs_native, ...
    cfg.array.center_freq, cfg.array.n_cycles, 'Plot', false));

% ---- 传感器：32 阵元同时接收 ----
sensor = struct();
sensor.mask = false(Nx, Ny);
sensor.mask(S.array_depth_idx, S.array_y_idx) = true;
sensor.record = {'p'};
n_rx = nnz(sensor.mask);
if n_rx ~= n_tx
    error('skull_simulation_function:sensor_mismatch', ...
        '阵元在网格上发生重叠：%d 个阵元仅占 %d 个网格点。', n_tx, n_rx);
end

% ---- 配置记录 ----
% 注：阵元位置、网格、时间步等全部在 sim_config 里，保存出的变量集合与
%     simdata_32el_15MHz_2um_curved_skull_sphere_<before|after>.mat 完全一致。
sim_config = struct();
sim_config.case_tag = case_def.tag;
sim_config.case_short = case_def.short;
sim_config.case_label = case_def.label;
sim_config.c_gel = cfg.media.c_gel;
sim_config.rho_gel = cfg.media.rho_gel;
sim_config.c_skull = case_def.c_skull;
sim_config.rho_skull = case_def.rho_skull;
sim_config.c_brain = cfg.media.c_brain;
sim_config.rho_brain = cfg.media.rho_brain;
sim_config.gel_thickness = cfg.geometry.gel_thickness;
sim_config.skull_thickness = cfg.geometry.skull_thickness;
sim_config.brain_thickness = cfg.geometry.brain_thickness;
sim_config.radius_gel_skull = cfg.geometry.radius_gel_skull;
sim_config.radius_skull_brain = cfg.geometry.radius_skull_brain;
sim_config.center_freq = cfg.array.center_freq;
sim_config.n_cycles = cfg.array.n_cycles;
sim_config.source_pressure = cfg.array.source_pressure;
sim_config.n_elements = n_tx;
sim_config.pitch = cfg.array.pitch;
sim_config.array_depth = S.array_depth;
sim_config.array_y_requested = S.array_y_requested;
sim_config.array_y_actual = S.array_y_actual;
sim_config.defects = cfg.defects;
sim_config.n_defects = numel(geom.defects);
sim_config.defect_geom = rmfield(geom.defects, {'label'});
sim_config.dx = cfg.grid.dx;
sim_config.dy = S.dy;
sim_config.dt = dt;                       % 仿真时间步进（不变）
sim_config.fs_native = fs_native;         % 仿真内部采样率 1/dt ≈ 4.47 GHz
sim_config.fs_out = fs_out;               % 导出数据采样率（300 MHz）
sim_config.dt_out = 1 / fs_out;
sim_config.Nt_native = kgrid.Nt;          % 仿真步数（动画用它）
sim_config.t_end = t_axis(end);           % 导出数据时长（8 us）
sim_config.Nt = n_out;                    % 导出采样点数
sim_config.Nx = Nx;
sim_config.Ny = Ny;
sim_config.domain_depth = Nx * cfg.grid.dx;
sim_config.lateral_width = Ny * S.dy;
sim_config.pml_size = cfg.grid.pml_size;
sim_config.pml_alpha = cfg.grid.pml_alpha;
sim_config.rf_size = [n_out, n_rx, n_tx];
sim_config.rf_dim_order = 'rf(t, rx, tx)';
sim_config.rf_units = 'Pa';
sim_config.z_gel_bottom = geom.z_gel_bottom;
sim_config.z_skull_bottom = geom.z_skull_bottom;
sim_config.z_brain_bottom = geom.z_brain_bottom;
sim_config.mat_file = mat_file;
sim_config.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

fprintf('  仿真步进 dt = %.4g s（%.3g Hz，%d 步）；导出采样 fs = %.4g Hz（%d 点）\n', ...
    dt, fs_native, kgrid.Nt, fs_out, n_out);
fprintf('  DataCast = %s\n', data_cast);

% ---- 断点续算 ----
ckpt_file = strrep(mat_file, '_full_data.mat', '_full_data_checkpoint.mat');
rf = [];
tx_start = 1;
if cfg.sim.resume && exist(ckpt_file, 'file')
    try
        L = load(ckpt_file, 'rf', 'tx_done');
        if isequal(size(L.rf), [kgrid.Nt, n_rx, n_tx]) && L.tx_done >= 1
            rf = L.rf;
            tx_start = L.tx_done + 1;
            fprintf('  从断点继续：已完成 %d/%d 个阵元\n', L.tx_done, n_tx);
        else
            fprintf('  断点文件与当前网格不一致，重新开始\n');
        end
    catch ME
        warning('skull_simulation_function:checkpoint', ...
            '读取断点失败（%s），重新开始', ME.message);
    end
end
if isempty(rf)
    rf = zeros(kgrid.Nt, n_rx, n_tx, 'single');
end
if tx_start > n_tx
    fprintf('  全部 %d 个阵元已完成，跳过采集\n', n_tx);
end

% ---- 逐阵元发射 ----
t_elapsed = 0;
for tx = tx_start:n_tx
    t_tx = tic;

    source = struct();
    source.p_mask = false(Nx, Ny);
    source.p_mask(S.array_depth_idx, S.array_y_idx(tx)) = true;
    source.p = input_signal;

    opts = kwave_opts(cfg, data_cast);
    try
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, opts{:});
    catch ME
        if strcmp(data_cast, 'gpuArray-single')
            warning('skull_simulation_function:gpu', ...
                'GPU 运行失败（%s），改用 CPU 重试。', ME.message);
            data_cast = 'single';
            opts = kwave_opts(cfg, data_cast);
            sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, opts{:});
        else
            rethrow(ME);
        end
    end

    if strcmp(data_cast, 'gpuArray-single')
        p_rx = gather(sensor_data.p);
    else
        p_rx = sensor_data.p;
    end
    rf(:, :, tx) = single(reshape(p_rx, n_rx, kgrid.Nt)).';

    el = toc(t_tx);
    t_elapsed = t_elapsed + el;
    if cfg.sim.verbose
        if tx < n_tx
            fprintf('  发射 %02d/%02d 完成，本阵元 %.1f s，预计剩余 %.1f min\n', ...
                tx, n_tx, el, (n_tx - tx) * (t_elapsed / max(1, tx - tx_start + 1)) / 60);
        else
            fprintf('  发射 %02d/%02d 完成，本阵元 %.1f s\n', tx, n_tx, el);
        end
    end

    % 断点：每 N 个阵元写一次
    if cfg.sim.save_checkpoint && tx < n_tx && ...
            mod(tx, cfg.sim.checkpoint_every) == 0
        tx_done = tx; %#ok<NASGU>
        save(ckpt_file, 'rf', 'tx_done', '-v7.3');
        fprintf('  已写入断点文件（%d/%d）\n', tx, n_tx);
    end
end

sim_config.data_cast = data_cast;
sim_config.elapsed = toc(t_start);

% ---- 重采样到导出采样率 fs_out（仿真步进 dt 不变，只是保存前重取样）----
if n_out ~= kgrid.Nt
    t_rs = tic;
    rf_out = interp1(t_native(:), reshape(rf, [kgrid.Nt, n_rx * n_tx]), t_axis, 'linear');
    rf = reshape(single(rf_out), [n_out, n_rx, n_tx]);
    clear rf_out
    fprintf('  已重采样：%d 点（%.3g Hz）-> %d 点（%.3g Hz），用时 %.1f s\n', ...
        kgrid.Nt, fs_native, n_out, fs_out, toc(t_rs));
end

% ---- 保存相控阵数据（3_signal，变量与 simdata_*.mat 完全一致）----
if cfg.sim.save_full_data
    save(mat_file, 'rf', 'sim_config', 'medium', 'input_signal', ...
        't_axis', 'x_axis', 'y_axis', '-v7.3');
    if exist(ckpt_file, 'file')
        delete(ckpt_file);
    end
    finfo = dir(mat_file);
    fprintf('  已保存相控阵数据：%s（%.1f MB）\n', mat_file, finfo.bytes / 2^20);
end

% 动画仍按仿真步进推进，所以返回原生时间轴
t_axis = t_native; %#ok<NASGU>
end


%% =====================================================================
%% 声场传播动画（5_video）
%% =====================================================================
function movie_file = run_animation(cfg, case_def, S, medium, sim_config, input_signal, t_axis, video_dir)
%RUN_ANIMATION 录制声场传播动画
%
%   额外跑一次仿真：激励孔径由 cfg.animation.aperture 决定——
%     'center'（当前默认）只激励中心阵元，'full' 为 32 个阵元同时同相激励；
%   用 k-Wave 的 PlotSim + RecordMovie 录制 p 场随时间的演化。
%
%   注意：
%     * 录制需要显示图窗，请在 MATLAB 桌面环境运行（不要用 -batch 无界面方式）；
%     * 帧数 = Nt / plot_freq，plot_freq 越小动画越平滑、录制越慢、文件越大；
%     * 坐标格式（x / mm、z / mm、深度从 0 起、Times New Roman）来自本地
%       k-Wave 副本里 kspaceFirstOrder2D.m 的 "[本工程修改]" 段落。

Nx = S.Nx;
Ny = S.Ny;
data_cast = S.data_cast;
n_el = cfg.array.n_elements;

movie_base = fullfile(video_dir, sprintf('%s_propagation_movie', case_def.short));

% ---- 激励孔径 ----
source = struct();
source.p_mask = false(Nx, Ny);
switch lower(cfg.animation.aperture)
    case 'center'
        ic = round(n_el / 2);
        source.p_mask(S.array_depth_idx, S.array_y_idx(ic)) = true;
        source.p = input_signal;
        aperture_txt = sprintf('中心阵元 %d', ic);
    otherwise   % 'full'
        source.p_mask(S.array_depth_idx, S.array_y_idx) = true;
        source.p = repmat(input_signal(:).', n_el, 1);
        aperture_txt = '全孔径 32 阵元同时激励';
end

% ---- 传感器：只取一个点，仅为了满足 k-Wave 的输入要求 ----
sensor = struct();
sensor.mask = false(Nx, Ny);
sensor.mask(S.array_depth_idx, S.array_y_idx(round(n_el / 2))) = true;
sensor.record = {'p'};

kgrid = kWaveGrid(Nx, cfg.grid.dx, Ny, S.dy);
kgrid.dt = sim_config.dt;
kgrid.Nt = numel(t_axis);

extra = {'PlotSim', true, ...
    'PlotFreq', cfg.animation.plot_freq, ...
    'PlotScale', cfg.animation.plot_scale, ...
    'RecordMovie', true, ...
    'MovieName', movie_base, ...
    'MovieProfile', cfg.animation.profile, ...
    'MovieArgs', {'FrameRate', cfg.animation.frame_rate}};
opts = kwave_opts(cfg, data_cast, extra);

fprintf('  录制声场动画：%s（%s，每 %d 步一帧，共约 %d 帧）\n', ...
    [movie_base '.avi'], aperture_txt, cfg.animation.plot_freq, ...
    ceil(numel(t_axis) / cfg.animation.plot_freq));

% k-Wave 内部用 figure() 建动画图窗，这里临时改默认图窗尺寸，
% 让竖长视场（3 mm x 6.25 mm）占据更大画面；函数退出时自动恢复。
old_fig_pos = get(groot, 'defaultFigurePosition');
set(groot, 'defaultFigurePosition', cfg.animation.figure_position);
restore_pos = onCleanup(@() set(groot, 'defaultFigurePosition', old_fig_pos)); %#ok<NASGU>

t_movie = tic;
try
    kspaceFirstOrder2D(kgrid, medium, source, sensor, opts{:});
catch ME
    warning('skull_simulation_function:movie', ...
        '动画录制失败（%s）。若是无界面运行请改为在 MATLAB 桌面运行。', ME.message);
    movie_file = '';
    return
end

cand = {[movie_base '.avi'], movie_base};
movie_file = '';
for i = 1:numel(cand)
    if exist(cand{i}, 'file')
        movie_file = cand{i};
        break;
    end
end
if isempty(movie_file)
    warning('skull_simulation_function:movie', '未找到生成的动画文件：%s', movie_base);
else
    finfo = dir(movie_file);
    fprintf('  已保存声场动画：%s（%.1f MB，用时 %.1f min）\n', ...
        movie_file, finfo.bytes / 2^20, toc(t_movie) / 60);
end
end

