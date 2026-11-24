function out = skull_signal_postprocess(cfg, opts)
%SKULL_SIGNAL_POSTPROCESS 颅骨仿真数据后处理：单阵元信号 + DAS（TFM）成像 + 中心纵向幅值分布
%
%  文件 3/3：skull_signal_postprocess.m
%
%  输入数据（按优先级自动查找，只需其中一套）：
%    1) 3_signal\<tag>_full_data.mat      仿真阶段（precheck = 1）保存的相控阵数据
%       （变量与 simdata_32el_15MHz_2um_curved_skull_sphere_<before|after>.mat 一致）
%    2) simdata_32el_15MHz_2um_curved_skull_sphere_<before|after>.mat
%       （根目录 / old / 1_code 下都会找）
%    3) <before|after>.mat               只有 rf 变量的旧数据（几何取配置）
%    tag = original（OACA 处理前）/ OACA（OACA 处理后）
%
%  输出（文件夹 3_signal 与 4_Image，命名严格按规范）：
%    3_signal\  OACA_signal.fig / .png       中心单阵元回波信号（处理前后同一张图：
%                                            处理前黑色虚线、处理后红色实线）
%               <tag>_signal_data.mat         波形数据（t、waveform、envelope...）
%               <tag>_full_data.mat           相控阵检测数据（rf/sim_config/medium/...）
%    4_Image\    <tag>_image.fig / .png       DAS 成像（dB）
%               <tag>_image_data.mat          图像数据（x、z、image、envelope...）
%               <tag>_axial_curve.fig / .png  沿中心的纵向幅值分布
%               <tag>_axial_curve.mat         曲线数据
%
%  用法：
%    skull_signal_postprocess                 % 自动读取配置与现有数据
%    skull_signal_postprocess(CFG)            % 指定配置
%    skull_signal_postprocess(CFG, opts)      % 附加控制选项
%
%    opts.cases                   需要处理的 tag 元胞数组（默认 cfg.run_cases）
%    opts.copy_existing_full_data 旧数据是否复制为 <tag>_full_data.mat（默认 true）
%
%  技术要点：
%    * 相控阵检测参数与 skull_simulation_model.m 第 5 节共用；
%    * 成像介质声速取大脑声速（cfg.media.c_brain）；
%    * 采样率取仿真数据的时间步进 fs = 1/dt；
%    * 成像网格由 cfg.imaging 的 x_range / z_range / step 手动设置；
%    * 成像深度范围默认取到 6.0 mm（计算域 6.25 mm 以内），
%      再往下是 PML 吸收层，延时叠加会在那里形成边界伪影，不是真实回波；
%    * 核心 TFM 成像模块为本文件末尾的 tfm_reconstruct_cpu（纯 CPU、向量化）。
%
%  参考：E:\3_HIS_APP_V1.0\+his\+algo\+tfm（TimeOfFlight.m / Interp.m / Algorithm.m）
% =========================================================================

%% ---------------- 0. 输入解析与配置 ----------------
if nargin < 1
    cfg = [];
end
if nargin < 2 || isempty(opts)
    opts = struct();
end

cfg = resolve_cfg(cfg);
cfg = local_normalize(cfg);

if ~isfield(opts, 'cases') || isempty(opts.cases)
    opts.cases = cfg.run_cases;
end
if ~isfield(opts, 'copy_existing_full_data') || isempty(opts.copy_existing_full_data)
    opts.copy_existing_full_data = true;
end

%% ---------------- 1. 输出目录 ----------------
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
out = struct();
out.folder = folder;

%% ---------------- 2. 成像网格 ----------------
X = cfg.imaging.x_range(1):cfg.imaging.step:cfg.imaging.x_range(2);
Z = cfg.imaging.z_range(1):cfg.imaging.step:cfg.imaging.z_range(2);

fprintf('\n================ 后处理：单阵元信号 + DAS 成像 ================\n');
fprintf('成像网格：x %s mm（%d 点），z %s mm（%d 点），步长 %.3f mm\n', ...
    mat2str([X(1) X(end)] * 1e3), numel(X), mat2str([Z(1) Z(end)] * 1e3), ...
    numel(Z), cfg.imaging.step * 1e3);
fprintf('成像声速：%.0f m/s（大脑声速），显示动态范围 %d dB\n', ...
    cfg.media.c_brain, cfg.imaging.db_range);
fprintf('==============================================================\n');

%% ---------------- 3. 逐状态处理 ----------------
metrics = struct('tag', {}, 'label', {}, 'echo_time', {}, 'echo_peak', {}, ...
    'early_peak', {}, 'defect_peak_depth', {}, 'defect_peak_db', {}, ...
    'axial_fwhm', {}, 'lateral_fwhm', {}, 'source_file', {});
file_list = struct('tag', {}, 'files', {});
sigs = [];                       % 收集各状态的单阵元信号，最后画一张对比图

for ic = 1:numel(opts.cases)
    tag = opts.cases{ic};
    k = find(strcmp({cfg.cases.tag}, tag), 1);
    if isempty(k)
        warning('skull_signal_postprocess:unknown_case', ...
            '未定义的配置标签 %s，已跳过。', tag);
        continue;
    end
    case_def = cfg.cases(k);

    fprintf('\n---------------- [%d/%d] %s ----------------\n', ...
        ic, numel(opts.cases), case_def.label);

    % ---- 3.1 读取相控阵数据 ----
    [D, used_file] = load_case_data(cfg, case_def, folder, opts);

    % ---- 3.2 中心单阵元信号数据（3_signal\<tag>_signal_data.mat）----
    sig = process_signal_data(cfg, case_def, D, folder);
    if isempty(sigs)
        sigs = sig;
    else
        sigs(end + 1) = sig; %#ok<SAGROW>
    end

    % ---- 3.3 DAS 成像 + 中心纵向幅值分布（4_Image）----
    img = process_image(cfg, case_def, D, folder, X, Z);

    % ---- 3.4 记录指标 ----
    metrics(end + 1) = struct('tag', case_def.tag, 'label', case_def.label, ...
        'echo_time', sig.echo_peak_time, 'echo_peak', sig.echo_peak_value, ...
        'early_peak', sig.early_peak_value, ...
        'defect_peak_depth', img.axial.defect_peak_depth, ...
        'defect_peak_db', img.axial.defect_peak_db, ...
        'axial_fwhm', img.axial.fwhm, ...
        'lateral_fwhm', img.lateral_fwhm, 'source_file', used_file); %#ok<SAGROW>
    file_list(end + 1) = struct('tag', case_def.tag, ...
        'files', {[{sig.mat_file}, img.files]}); %#ok<SAGROW>
end

% ---- 3.5 处理前后信号对比图（一张图，放在 3_signal）----
if ~isempty(sigs)
    sig_paths = plot_signal_figure(cfg, sigs, folder);
    file_list(end + 1) = struct('tag', 'signal-figure', 'files', {sig_paths});
end

out.metrics = metrics;
out.files = file_list;

%% ---------------- 4. 汇总 ----------------
fprintf('\n======================== 后处理汇总 ========================\n');
fprintf('%-14s %9s %10s %10s %10s %10s\n', ...
    '处理状态', '回波/us', '回波-近场', '缺陷峰/mm', '缺陷峰/dB', '纵向-6dB/mm');
for i = 1:numel(metrics)
    m = metrics(i);
    fprintf('%-14s %9.3f %9.1f dB %9.3f %9.1f %9.3f\n', ...
        m.label, m.echo_time * 1e6, ...
        20 * log10(max(m.echo_peak, eps) / max(m.early_peak, eps)), ...
        m.defect_peak_depth * 1e3, m.defect_peak_db, m.axial_fwhm * 1e3);
end
fprintf('横向 -6 dB 宽度：');
for i = 1:numel(metrics)
    fprintf('%s %.3f mm  ', metrics(i).tag, metrics(i).lateral_fwhm * 1e3);
end
fprintf('\n');
fprintf('----------------------------------------------------\n');
for i = 1:numel(file_list)
    fprintf('[%s] 共生成 %d 个文件\n', file_list(i).tag, numel(file_list(i).files));
    for j = 1:numel(file_list(i).files)
        fprintf('    %s\n', file_list(i).files{j});
    end
end
fprintf('====================================================\n');

end % ==================== 主函数结束 ====================


%% =====================================================================
%% 配置读取
%% =====================================================================
function cfg = resolve_cfg(cfg)
%RESOLVE_CFG 配置读取顺序：入参 -> 基础工作区 CFG -> 配置快照 mat -> 内置默认值

if ~isempty(cfg)
    return
end

% 1) 基础工作区里的 CFG（运行过 skull_simulation_model.m 就会有）
try
    if evalin('base', 'exist(''CFG'',''var'')') == 1
        cfg = evalin('base', 'CFG');
        fprintf('已读取基础工作区中的配置 CFG\n');
        return
    end
catch
    % 忽略：无基础工作区访问权限时继续往下找
end

% 2) 配置快照
snap = fullfile(fileparts(mfilename('fullpath')), 'skull_model_cfg.mat');
if exist(snap, 'file')
    L = load(snap, 'CFG');
    if isfield(L, 'CFG')
        cfg = L.CFG;
        fprintf('已读取配置快照：%s\n', snap);
        return
    end
end

% 3) 内置默认值（与 skull_simulation_model.m 保持一致）
cfg = default_cfg();
warning('skull_signal_postprocess:default_cfg', ...
    ['未找到配置（skull_model_cfg.mat），已使用内置默认值。\n' ...
    '建议先运行 skull_simulation_model.m（precheck = 0 即可）生成配置快照。']);
end


function cfg = default_cfg()
%DEFAULT_CFG 内置默认配置（与 skull_simulation_model.m 保持一致）
here = fileparts(mfilename('fullpath'));          % 代码目录（1_code）
cfg = struct();
cfg.precheck = 0;
cfg.run_cases = {'original', 'OACA'};
cfg.code_dir = here;
cfg.root_dir = fileparts(here);                   % 工程根目录 = 代码目录的上一级
if isempty(cfg.root_dir) || ~exist(cfg.root_dir, 'dir')
    cfg.root_dir = here;
end
cfg.kwave_dir = fullfile(here, 'k-wave-toolbox-version-1.4', 'k-Wave');
cfg.folder = struct('model', '2_skull_model', 'signal', '3_signal', ...
    'image', '4_Image', 'video', '5_video');
cfg.media = struct('c_gel', 1490, 'rho_gel', 998, ...
    'c_brain', 1540, 'rho_brain', 1050);
cfg.cases = struct( ...
    'tag', {'original', 'OACA'}, ...
    'label', {'OACA 处理前', 'OACA 处理后'}, ...
    'short', {'original', 'OACA'}, ...
    'c_skull', {3110, 1680}, ...
    'rho_skull', {1975, 1152});
cfg.array = struct('n_elements', 32, 'pitch', 50e-6, 'center_freq', 15e6, ...
    'n_cycles', 1, 'source_pressure', 1e3, 'element_width', 50e-6);
cfg.grid = struct('dx', 1e-6, 'lateral_width', 3.002e-3, 'dt_div', 3, ...
    't_end', 8e-6, 'fs_out', 300e6, 'pml_size', 50, 'pml_alpha', 1e3);
cfg.imaging = struct('x_range', [-1.5e-3, 1.5e-3], 'z_range', [0, 6.0e-3], ...
    'step', 0.02e-3, 'db_range', 40, 'colormap', 'jet', ...
    'center_element', [], 'axial_lateral', 0);
cfg.output = struct('model_figures', true, 'export_png', true, ...
    'close_figures', true);
cfg.model_plot = struct('clim_speed', [1400 3400], 'clim_density', [900 2000], ...
    'font_name', 'Times New Roman', 'font_size', 12, ...
    'figure_size', [277.6 404], 'draw_impedance', false);
cfg.defects = struct('depth', 5.0e-3, 'lateral', 0, 'diameter', 0.1e-3, ...
    'c', 2000, 'rho', 1500, 'label', '散射体');
end


function cfg = local_normalize(cfg)
%LOCAL_NORMALIZE 补齐缺省字段

d_img = struct('x_range', [-1.5e-3, 1.5e-3], 'z_range', [0, 6.0e-3], ...
    'step', 0.02e-3, 'db_range', 40, 'colormap', 'jet', ...
    'center_element', [], 'axial_lateral', 0);
if ~isfield(cfg, 'imaging') || isempty(cfg.imaging)
    cfg.imaging = d_img;
else
    f = fieldnames(d_img);
    for i = 1:numel(f)
        if ~isfield(cfg.imaging, f{i}) || isempty(cfg.imaging.(f{i}))
            cfg.imaging.(f{i}) = d_img.(f{i});
        end
    end
end

d_out = struct('model_figures', true, 'export_png', true, ...
    'close_figures', true);
if ~isfield(cfg, 'output') || isempty(cfg.output)
    cfg.output = d_out;
else
    f = fieldnames(d_out);
    for i = 1:numel(f)
        if ~isfield(cfg.output, f{i}) || isempty(cfg.output.(f{i}))
            cfg.output.(f{i}) = d_out.(f{i});
        end
    end
end

if ~isfield(cfg, 'run_cases') || isempty(cfg.run_cases)
    cfg.run_cases = {cfg.cases.tag};
end
if ~isfield(cfg, 'defects') || isempty(cfg.defects)
    cfg.defects = struct('depth', 5.0e-3, 'lateral', 0, 'diameter', 0.1e-3, ...
        'c', 2000, 'rho', 1500, 'label', '散射体');
end
if ~isfield(cfg, 'folder') || isempty(cfg.folder)
    cfg.folder = struct('model', '2_skull_model', 'signal', '3_signal', ...
        'image', '4_Image', 'video', '5_video');
end
if ~isfield(cfg, 'code_dir') || isempty(cfg.code_dir)
    cfg.code_dir = fileparts(mfilename('fullpath'));
end
if ~isfield(cfg, 'root_dir') || isempty(cfg.root_dir)
    cfg.root_dir = fileparts(cfg.code_dir);
    if isempty(cfg.root_dir) || ~exist(cfg.root_dir, 'dir')
        cfg.root_dir = cfg.code_dir;
    end
end
end


%% =====================================================================
%% 数据读取
%% =====================================================================
function [D, used_file] = load_case_data(cfg, case_def, folder, opts)
%LOAD_CASE_DATA 读取某个处理状态的相控阵数据，并整理成统一结构
%
%   优先读取 3_signal\<tag>_full_data.mat；没有则回退到工程根目录下已有的
%   simdata_*.mat（或 <before|after>.mat），便于先用现有数据检查后处理效果。

target = fullfile(folder.signal, sprintf('%s_full_data.mat', case_def.short));
cands = candidate_files(cfg, case_def);

used_file = '';
for i = 1:numel(cands)
    if exist(cands{i}, 'file')
        used_file = cands{i};
        break;
    end
end
if isempty(used_file)
    error('skull_signal_postprocess:no_data', ...
        ['未找到 %s 的相控阵数据。已尝试：\n%s\n' ...
        '请运行 skull_simulation_function（precheck = 1）生成数据。'], ...
        case_def.label, strjoin(cands, newline));
end

fprintf('  数据来源：%s\n', used_file);
t0 = tic;

% 只读需要的变量，避免把 140 MB 的 medium 也读进来
w = whos('-file', used_file);
have = {w.name};
want = intersect({'rf', 'sim_config', 't_axis', 'x_axis', 'y_axis', ...
    'array_info', 'input_signal'}, have);
S = load(used_file, want{:});

D = parse_data(S, cfg, case_def);
D.source_file = used_file;
fprintf('  已读取 rf：%s（%s），用时 %.1f s\n', ...
    mat2str(size(D.rf)), class(D.rf), toc(t0));

% 旧数据 -> 复制为规范命名的相控阵数据文件
if opts.copy_existing_full_data && ~strcmpi(used_file, target) && ~exist(target, 'file')
    try
        copyfile(used_file, target);
        fprintf('  已复制为规范命名：%s\n', target);
    catch ME
        warning('skull_signal_postprocess:copy', ...
            '复制相控阵数据失败（%s），已跳过', ME.message);
    end
end
end


function cands = candidate_files(cfg, case_def)
%CANDIDATE_FILES 候选数据文件（按优先级）
cands = {fullfile(cfg.root_dir, cfg.folder.signal, ...
    sprintf('%s_full_data.mat', case_def.short))};

% 旧命名：simdata_<stem>_<before|after>.mat（可能放在根目录 / old / 代码目录）
legacy = legacy_tag(case_def.tag);
stem = 'simdata_32el_15MHz_2um_curved_skull_sphere';
roots = {cfg.root_dir, cfg.code_dir, ...
    fullfile(cfg.root_dir, 'old'), fullfile(cfg.code_dir, 'old')};
for i = 1:numel(roots)
    cands{end + 1} = fullfile(roots{i}, sprintf('%s_%s.mat', stem, legacy)); %#ok<AGROW>
    cands{end + 1} = fullfile(roots{i}, sprintf('%s.mat', legacy)); %#ok<AGROW>
end

% 兜底：工程根目录下递归找 simdata_*_<tag>*.mat
d = dir(fullfile(cfg.root_dir, '**', sprintf('simdata_*_%s*.mat', legacy)));
for i = 1:numel(d)
    cands{end + 1} = fullfile(d(i).folder, d(i).name); %#ok<AGROW>
end
end


function t = legacy_tag(tag)
%LEGACY_TAG 新命名 -> 旧命名
switch lower(tag)
    case 'original'
        t = 'before';
    case 'oaca'
        t = 'after';
    otherwise
        t = lower(tag);
end
end


function D = parse_data(S, cfg, case_def)
%PARSE_DATA 把不同来源的数据整理成统一结构

D = struct();
if ~isfield(S, 'rf')
    error('skull_signal_postprocess:no_rf', ...
        '数据文件中没有 rf 变量（相控阵数据）：%s', case_def.label);
end

rf = S.rf;
D.rf = rf;
[D.Nt, D.nRx, D.nTx] = size(rf);

sc = [];
if isfield(S, 'sim_config')
    sc = S.sim_config;
end

% ---- 时间轴与采样率 ----
if isfield(S, 't_axis') && numel(S.t_axis) == D.Nt
    D.t = double(S.t_axis(:).');
    D.fs = 1 / (D.t(2) - D.t(1));
else
    if ~isempty(sc) && isfield(sc, 'dt')
        D.fs = 1 / sc.dt;
    else
        error('skull_signal_postprocess:no_time', '数据中缺少时间轴 t_axis / dt。');
    end
    D.t = (0:D.Nt - 1) / D.fs;
end

% ---- 介质声速（成像介质 = 大脑）----
D.c_brain = cfg.media.c_brain;
if ~isempty(sc) && isfield(sc, 'c_brain')
    D.c_brain = sc.c_brain;
end

% ---- 阵元位置 ----
D.pitch = cfg.array.pitch;
if ~isempty(sc) && isfield(sc, 'pitch')
    D.pitch = sc.pitch;
end
if ~isempty(sc) && isfield(sc, 'array_y_actual')
    D.array_y = double(sc.array_y_actual(:).');
else
    D.array_y = ((1:D.nTx) - (D.nTx + 1) / 2) * D.pitch;
end
if ~isempty(sc) && isfield(sc, 'array_depth')
    D.array_depth = double(sc.array_depth);
else
    D.array_depth = cfg.grid.dx;
end
if isfield(S, 'array_info') && isfield(S.array_info, 'y_actual')
    D.array_y = double(S.array_info.y_actual(:).');
end

% ---- 缺陷位置 ----
D.defects = cfg.defects;
if ~isempty(sc) && isfield(sc, 'defects') && ~isempty(sc.defects)
    D.defects = sc.defects;
elseif ~isempty(sc) && isfield(sc, 'scatterer_depth')
    D.defects = struct('depth', sc.scatterer_depth, 'lateral', sc.scatterer_y, ...
        'diameter', sc.scatterer_diameter, 'c', sc.scatterer_c, ...
        'rho', sc.scatterer_rho, 'label', '散射体');
end

D.case_tag = case_def.tag;
D.case_short = case_def.short;
D.case_label = case_def.label;
end


%% =====================================================================
%% 中心单阵元信号波形（文件夹 3_signal）
%% =====================================================================
function sig = process_signal_data(cfg, case_def, D, folder)
%PROCESS_SIGNAL_DATA 提取中心阵元的自发自收信号，计算指标并保存 mat 数据
%
%   图统一在 plot_signal_figure 里画（处理前后放同一张图），这里只负责数据。

ic = cfg.imaging.center_element;
if isempty(ic)
    ic = round(D.nTx / 2);
end
ic = max(1, min(min(D.nTx, D.nRx), ic));

wave = double(D.rf(:, ic, ic));
wave = wave(:);
env = envelope_hilbert(wave);

t = D.t(:);

% 缺陷名义回波时刻（双程）
if isempty(D.defects)
    echo_t = [];
else
    echo_t = 2 * [D.defects.depth] / D.c_brain;
end
if isempty(echo_t)
    echo_t = t(end) / 2;
end
t_echo = min(echo_t);

fprintf('  单阵元：第 %d 号（y = %.1f um），fs = %.4g Hz\n', ...
    ic, D.array_y(ic) * 1e6, D.fs);
fprintf('  缺陷名义回波时刻：%.3f us\n', t_echo * 1e6);

% ---- 回波峰值：在名义回波时刻 ±0.6 us 内取包络最大值 ----
win = [max(t(1), t_echo - 0.6e-6), min(t(end), t_echo + 0.6e-6)];
sel = (t >= win(1)) & (t <= win(2));
if ~any(sel)
    sel = true(size(t));
end
idx_sel = find(sel);
[echo_peak_value, ipk_rel] = max(env(idx_sel));
ipk = idx_sel(ipk_rel);
echo_peak_time = t(ipk);

% 近场（回波窗之前）的包络最大值，用于评估信杂比
early = t < win(1);
if any(early)
    early_peak_value = max(env(early));
else
    early_peak_value = max(env);
end

% ---- 保存数据 ----
signal_data = struct();
signal_data.case_tag = case_def.tag;
signal_data.case_short = case_def.short;
signal_data.case_label = case_def.label;
signal_data.element = ic;
signal_data.element_y = D.array_y(ic);
signal_data.t = t;
signal_data.fs = D.fs;
signal_data.waveform = wave;
signal_data.envelope = env;
signal_data.envelope_db = 20 * log10(max(env, eps) / max(max(env), eps));
signal_data.units = 'Pa';
signal_data.medium_speed = D.c_brain;
signal_data.defects = D.defects;
signal_data.echo_times = echo_t(:).';
signal_data.echo_peak_time = echo_peak_time;
signal_data.echo_peak_value = echo_peak_value;
signal_data.early_peak_value = early_peak_value;
signal_data.source_file = D.source_file;
signal_data.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

mat_file = fullfile(folder.signal, sprintf('%s_signal_data.mat', case_def.short));
save(mat_file, 'signal_data', '-v7.3');
fprintf('  已保存：%s\n', mat_file);

sig = struct();
sig.tag = case_def.tag;
sig.label = case_def.label;
sig.element = ic;
sig.element_y = D.array_y(ic);
sig.t = t;
sig.wave = wave;
sig.env = env;
sig.echo_t = echo_t(:).';
sig.t_echo = t_echo;
sig.echo_peak_time = echo_peak_time;
sig.echo_peak_value = echo_peak_value;
sig.early_peak_value = early_peak_value;
sig.defects = D.defects;
sig.mat_file = mat_file;
end


%% =====================================================================
%% 处理前后单阵元回波信号对比图（一张图）
%% =====================================================================
function paths = plot_signal_figure(cfg, sigs, folder)
%PLOT_SIGNAL_FIGURE 处理前后中心单阵元回波信号画在同一张图里
%
%   版式按 3_signal\OACA_signal.fig 定稿模板：
%     单坐标区（不再上下分栏，图窗按内容收紧，不留空白），
%     Times New Roman 12，网格开，粗线（LW 2），
%     处理前 = 黑色虚线，处理后 = 红色实线，名义回波时刻画竖虚线，
%     左上角图例标明 before / after OACA（不再另加底部图注）。

if isempty(sigs)
    paths = {};
    return
end

t_echo = min([sigs.t_echo]);

% ---- 图窗：按模板收紧高度，去掉底部空白（模板里那行字已取消）----
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 596.8 210]);
ax = axes(fig, 'Position', [0.105 0.20 0.875 0.70]);
hold(ax, 'on');
grid(ax, 'on');

% ---- 曲线：先处理前（黑虚线），再处理后（红实线）----
h = gobjects(1, numel(sigs));
for i = 1:numel(sigs)
    s = sigs(i);
    if strcmpi(s.tag, 'original')
        style = {'--', [0 0 0], 'before OACA'};
    else
        style = {'-', [1 0 0], 'after OACA'};
    end
    h(i) = plot(ax, s.t * 1e6, s.wave, style{1}, 'Color', style{2}, ...
        'LineWidth', 2, 'DisplayName', style{3});
end

% ---- 只显示缺陷回波附近 ±1.2 us ----
lo = max(sigs(1).t(1), t_echo - 1.2e-6) * 1e6;
hi = min(sigs(1).t(end), t_echo + 1.2e-6) * 1e6;
xlim(ax, [lo hi]);

% 纵轴按局部数据自适应（留 15% 余量）
ymax = 0;
for i = 1:numel(sigs)
    sel = sigs(i).t >= lo * 1e-6 & sigs(i).t <= hi * 1e-6;
    ymax = max(ymax, max(abs(sigs(i).wave(sel))));
end
if ~(ymax > 0)
    ymax = 1;
end
ylim(ax, [-1.15 1.15] * ymax);

% ---- 名义回波时刻 ----
xline(ax, t_echo * 1e6, 'k--', 'LineWidth', 0.5);

% ---- 坐标区版式（按 3_signal\OACA_signal.fig 定稿模板）----
set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'TickDir', 'in', 'Box', 'on');
title(ax, 'Echo Signal of the center element', 'FontSize', 12);
xlabel(ax, 't / \mus', 'FontSize', 12);
ylabel(ax, 'p / Pa', 'FontSize', 12);

lg = legend(ax, h, 'Location', 'northwest', 'FontSize', 10, 'Box', 'on');
set(lg, 'FontName', 'Times New Roman');

base = fullfile(folder.signal, 'OACA_signal');
paths = save_figure(fig, base, cfg.output.export_png, cfg.output.close_figures, 300);
fprintf('  已保存：%s.fig（处理前后对比，单图）\n', base);
end


%% =====================================================================
%% DAS（TFM 延时叠加）成像 + 中心纵向幅值分布（4_Image）
%% =====================================================================
function img = process_image(cfg, case_def, D, folder, X, Z)
%PROCESS_IMAGE DAS 成像、出图、中心纵向幅值分布

c = D.c_brain;
tx_y = D.array_y;
tx_z = D.array_depth * ones(size(tx_y));
rx_y = tx_y;
rx_z = tx_z;

fprintf('  DAS 成像：%d 发 x %d 收，网格 %d x %d，c = %.0f m/s\n', ...
    size(D.rf, 3), size(D.rf, 2), numel(Z), numel(X), c);

t_tfm = tic;
RFimg = tfm_reconstruct_cpu(D.rf, D.fs, tx_y, tx_z, rx_y, rx_z, X, Z, c);
env = envelope_hilbert(RFimg);
peak = max(env(:));
if ~(peak > 0)
    warning('skull_signal_postprocess:zero_image', ...
        '成像结果全为零，请检查数据与声速设置。');
    peak = 1;
end
img_db = db_normalize(env, cfg.imaging.db_range);
fprintf('  成像耗时 %.1f s，图像最大值 %.3f Pa\n', toc(t_tfm), peak);

% ---- 中心纵向幅值分布 ----
[~, ix0] = min(abs(X - cfg.imaging.axial_lateral));
curve = env(:, ix0);
curve_db = db_normalize(curve, cfg.imaging.db_range);
[pk_val, ipk] = max(curve);                 % 全局峰值（通常在阵元表面附近）

% 缺陷回波峰值：在缺陷深度 ±0.5 mm 内取最大值（比全局峰值更有物理意义）
z_all = Z(:);
if ~isempty(D.defects)
    [~, idef_lat] = min(abs([D.defects.lateral] - cfg.imaging.axial_lateral));
    z_def = D.defects(idef_lat).depth;
    def_win = abs(z_all - z_def) <= 0.5e-3;
    if ~any(def_win)
        def_win = true(size(z_all));
    end
    idx_win = find(def_win);
    [def_val, idef_rel] = max(curve(idx_win));
    idef = idx_win(idef_rel);
else
    z_def = NaN;
    idef = ipk;
    def_val = pk_val;
end

axial = struct();
axial.x_position = X(ix0);
axial.z = z_all;
axial.values = curve(:);
axial.values_db = curve_db(:);
axial.peak_depth = Z(ipk);                  % 全局峰值（表面）
axial.peak_value = pk_val;
axial.peak_value_db = curve_db(ipk);
axial.defect_depth = z_def;                 % 缺陷名义深度
axial.defect_peak_depth = Z(idef);          % 缺陷区实测峰值深度
axial.defect_peak_value = def_val;
axial.defect_peak_db = curve_db(idef);
axial.fwhm = fwhm_axis(z_all, curve(:), 6, idef);   % 缺陷回波 -6 dB 宽度
axial.defects = D.defects;
axial.case_label = case_def.label;
axial.source_file = D.source_file;

% ---- 横向 FWHM（在缺陷峰值深度处）----
lateral_fwhm = fwhm_axis(X(:).', env(idef, :), 6);

% ---- 出图：DAS 图像（版式按 4_Image\original_image.fig 模板）----
when = 'after';
if strcmpi(case_def.tag, 'original')
    when = 'before';
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 279.2 402.4]);
ax = axes(fig);
imagesc(ax, X * 1e3, Z * 1e3, img_db);
colormap(ax, cfg.imaging.colormap);
axis(ax, 'image');
set(ax, 'CLim', [-cfg.imaging.db_range 0], 'FontName', 'Times New Roman', ...
    'FontSize', 12, 'TickDir', 'in', 'Box', 'on');
xlabel(ax, 'x / mm');
ylabel(ax, 'z / mm');
title(ax, sprintf('DAS image %s NCADA', when), 'FontSize', 12);

cb = colorbar(ax);
set(cb, 'FontName', 'Times New Roman');
% 色标单位：图窗文本框，放在色标正下方（与模板一致）
cpos = get(cb, 'Position');
annotation(fig, 'textbox', [cpos(1), cpos(2) - 0.078, 0.132, 0.065], ...
    'String', 'dB', 'FontName', 'Times New Roman', 'FontSize', 10, ...
    'Color', [0 0 0], 'EdgeColor', 'none', 'BackgroundColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'Interpreter', 'tex', 'FitBoxToText', 'off');

base_img = fullfile(folder.image, sprintf('%s_image', case_def.short));
paths_img = save_figure(fig, base_img, cfg.output.export_png, cfg.output.close_figures, 300);

% ---- 出图：中心纵向幅值分布（版式按 4_Image\original_axial_curve.fig 模板）----
fig2 = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 413.6 191.2]);
ax2 = axes(fig2);
plot(ax2, Z * 1e3, curve_db, '-', 'Color', [0 0 0], 'LineWidth', 2);
grid(ax2, 'on');
hold(ax2, 'on');
xlim(ax2, [Z(1) Z(end)] * 1e3);
ylim(ax2, [-cfg.imaging.db_range 2]);
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'TickDir', 'in', 'Box', 'on');
xlabel(ax2, 'z / mm');
ylabel(ax2, 'Amplitude / dB');
title(ax2, sprintf('Axial amplitude curve at x=%g', X(ix0) * 1e3), 'FontSize', 12);

base_ax = fullfile(folder.image, sprintf('%s_axial_curve', case_def.short));
paths_ax = save_figure(fig2, base_ax, cfg.output.export_png, cfg.output.close_figures, 300);

% ---- 保存图像数据 ----
image_data = struct();
image_data.case_tag = case_def.tag;
image_data.case_short = case_def.short;
image_data.case_label = case_def.label;
image_data.x = X;
image_data.z = Z;
image_data.image = RFimg;              % 线性幅值（DAS 延时叠加结果）
image_data.envelope = env;             % 包络
image_data.image_db = img_db;          % 归一化 dB（下限 = -db_range）
image_data.db_range = cfg.imaging.db_range;
image_data.speed = c;
image_data.fs = D.fs;
image_data.n_elements = numel(tx_y);
image_data.pitch = D.pitch;
image_data.array_y = tx_y;
image_data.array_depth = D.array_depth;
image_data.nTx = size(D.rf, 3);
image_data.nRx = size(D.rf, 2);
image_data.Nt = size(D.rf, 1);
image_data.defects = D.defects;
image_data.axial_curve = axial;
image_data.lateral_fwhm = lateral_fwhm;
image_data.source_file = D.source_file;
image_data.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

mat_img = fullfile(folder.image, sprintf('%s_image_data.mat', case_def.short));
save(mat_img, 'image_data', '-v7.3');
mat_ax = fullfile(folder.image, sprintf('%s_axial_curve.mat', case_def.short));
save(mat_ax, 'axial', '-v7.3');

fprintf('  已保存：%s\n', [base_img '.fig']);
fprintf('  已保存：%s\n', mat_img);
fprintf('  已保存：%s\n', [base_ax '.fig']);
fprintf('  已保存：%s\n', mat_ax);

img = struct();
img.axial = axial;
img.lateral_fwhm = lateral_fwhm;
img.files = [paths_img, paths_ax, {mat_img}, {mat_ax}];
end


%% =====================================================================
%% 核心：TFM 全聚焦成像（CPU）
%% =====================================================================
function img = tfm_reconstruct_cpu(rf, fs, tx_y, tx_z, rx_y, rx_z, X, Z, c)
%TFM_RECONSTRUCT_CPU 全聚焦法（FMC）延时叠加成像，纯 CPU 实现
%
%   img(iz, ix) = sum_{tx, rx} rf( t = |p-tx|/c + |p-rx|/c, rx, tx )
%
%   输入：
%     rf     [Nt, nRx, nTx] 射频数据（Pa）
%     fs     采样率 [Hz]（= 1/dt，取仿真程序的时间步进）
%     tx_y   [1, nTx] 发射阵元横向坐标 [m]
%     tx_z   [1, nTx] 发射阵元深度坐标 [m]
%     rx_y   [1, nRx] 接收阵元横向坐标 [m]
%     rx_z   [1, nRx] 接收阵元深度坐标 [m]
%     X      [1, nX]  成像横向坐标 [m]
%     Z      [1, nZ]  成像深度坐标 [m]
%     c      成像介质声速 [m/s]（本工程取大脑声速）
%
%   输出：
%     img    [nZ, nX] 延时叠加结果（未取包络）
%
%   取数规则与 +his\+algo\+tfm 的 GPU 内核 / Interp.linear 一致：
%     sampleIndex = t * fs        （0 基）
%     仅当 0 <= sampleIndex < Nt-1 时参与叠加，越界按 0 处理
%     v = v0 + frac * (v1 - v0)   （线性插值）
%
%   实现要点：
%     * 声时表按"每个像素相对每个阵元"预先算好（发射/接收各一张），
%       通道循环里只做插值累加，避免重复开方；
%     * 全部 double 计算，结果与 GPU 内核逐点一致到单精度。

[Nt, nRx, nTx] = size(rf);
nX = numel(X);
nZ = numel(Z);
nPix = nX * nZ;

[Xg, Zg] = meshgrid(X, Z);          % Xg/Zg 均为 [nZ, nX]
px = Xg(:);
pz = Zg(:);

% ---- 声时表 ----
txT = zeros(nPix, nTx);
for it = 1:nTx
    txT(:, it) = sqrt((px - tx_y(it)).^2 + (pz - tx_z(it)).^2) / c;
end
rxT = zeros(nPix, nRx);
for ir = 1:nRx
    rxT(:, ir) = sqrt((px - rx_y(ir)).^2 + (pz - rx_z(ir)).^2) / c;
end

% ---- 延时叠加 ----
acc = zeros(nPix, 1);
last = Nt - 1;
for it = 1:nTx
    for ir = 1:nRx
        sample_time = txT(:, it) + rxT(:, ir);
        idx = sample_time * fs;                 % 0 基采样序号
        valid = (idx >= 0) & (idx < last);
        if ~any(valid)
            continue
        end
        idx_v = idx(valid);
        i0 = floor(idx_v) + 1;                  % 1 基
        frac = idx_v - floor(idx_v);
        trace = double(rf(:, ir, it));
        v0 = trace(i0);
        v1 = trace(i0 + 1);
        acc(valid) = acc(valid) + (v0 + frac .* (v1 - v0));
    end
end

img = reshape(acc, nZ, nX);
end


%% =====================================================================
%% 通用工具
%% =====================================================================
function env = envelope_hilbert(x, pad)
%ENVELOPE_HILBERT 沿第 1 维做 Hilbert 变换取包络（不依赖 Signal Processing Toolbox）
%
%   两端补零后再变换，抑制 FFT 循环卷积在数据边界产生的伪影：
%   直接对整列做 Hilbert 时，列的尾部会与首部"接上"，成像图底部的弱信号
%   会被放大成虚假的扇形强回波（实测可以放大 3 个数量级）。

n = size(x, 1);
if n < 2
    env = abs(x);
    return
end
if nargin < 2 || isempty(pad)
    pad = max(50, round(0.25 * n));
end

xp = [zeros(pad, size(x, 2), 'like', x); x; zeros(pad, size(x, 2), 'like', x)];
m = size(xp, 1);

X = fft(xp, [], 1);
h = zeros(m, 1);
if mod(m, 2) == 0
    h(1) = 1;
    h(m / 2 + 1) = 1;
    h(2:m / 2) = 2;
else
    h(1) = 1;
    h(2:(m + 1) / 2) = 2;
end
env = abs(ifft(X .* h, [], 1));
env = env(pad + 1:pad + n, :);
end


function db = db_normalize(v, db_range)
%DBNORMALIZE 归一化 dB（0 dB = 最大值），下限 -db_range
peak = max(v(:));
if ~(peak > 0) || ~isfinite(peak)
    db = -db_range * ones(size(v));
    return
end
db = 20 * log10(max(v, eps) / peak);
db = max(db, -db_range);
end


function w = fwhm_axis(ax, v, drop_db, ipk)
%FWHM_AXIS 计算 -drop_db dB 处的半高全宽（线性插值），返回物理宽度
%
%   ipk 可选：指定峰值所在下标（默认取全局最大值）
w = NaN;
v = v(:);
ax = ax(:);
if nargin < 4 || isempty(ipk)
    [pk, ipk] = max(v);
else
    pk = v(ipk);
end
if ~(pk > 0)
    return
end
thr = pk * 10^(-drop_db / 20);

% 左边界
i1 = find(v(1:ipk) < thr, 1, 'last');
if isempty(i1) || i1 >= numel(v)
    lo = ax(1);
else
    lo = interp1(v(i1:i1 + 1), ax(i1:i1 + 1), thr);
end

% 右边界
i2 = ipk - 1 + find(v(ipk:end) < thr, 1, 'first');
if isempty(i2) || i2 > numel(v) || i2 < 2
    hi = ax(end);
else
    hi = interp1(v(i2 - 1:i2), ax(i2 - 1:i2), thr);
end

if isfinite(lo) && isfinite(hi)
    w = abs(hi - lo);
end
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
        warning('skull_signal_postprocess:png', ...
            '导出 png 失败（%s），已跳过：%s.png', ME.message, base);
    end
end
if close_after
    close(fig);
end
end

