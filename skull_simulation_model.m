%% ==================== 颅骨超声相控阵仿真：模型配置 + 总控 ====================
%
%  文件 1/3：skull_simulation_model.m
%
%  作用：
%    1) 设置 OACA 处理前 / 处理后两套介质的声学参数（声速、密度）与几何参数；
%    2) 设置相控阵配置（阵元数、间距、中心频率、激励周期、发射声压）；
%    3) 设置缺陷（散射体）分布；
%    4) 设置是否保存声场传播动画，以及后处理所需的成像参数；
%    5) 作为整个工程的总入口：precheck = 0 时仅导出模型分布图，
%       precheck = 1 时跑完全仿真并接着做后处理。
%
%  运行方式：直接"运行"本文件（F5）。
%
%  ---------------------------------------------------------------------
%  precheck 按钮（见下面第 1 节）
%  ---------------------------------------------------------------------
%    precheck = 0 ：仅查看。只建立介质模型并导出 4 张 MATLAB Fig 图
%                   （2_skull_model：处理前后的声速、密度分布），不跑 k-Wave 仿真。
%                   用于快速确认几何、材料、缺陷摆放是否正确。
%    precheck = 1 ：跑完全仿真。在 precheck = 0 的基础上继续调用
%                     skull_simulation_function.m
%                       -> 逐阵元发射 / 全阵接收的 k-Wave 仿真
%                       -> 保存 OACA 处理前后的相控阵数据 mat 文件
%                       -> 保存 OACA 处理前后的声场传播动画 avi 文件
%                     skull_signal_postprocess.m
%                       -> 中心单阵元信号、DAS（TFM 延时叠加）成像、中心纵向幅值分布
%
%  ---------------------------------------------------------------------
%  输出目录（相对工程根目录 = 本文件所在目录 1_code 的上一级）
%  ---------------------------------------------------------------------
%    2_skull_model\ original_Speed       original_Density
%                   OACA_Speed           OACA_Density
%                   （各含 .fig，另附同名 .png 预览；版式按 OACA_Speed.fig 模板，
%                     处理前后色标范围统一：声速 1400~3400，密度 900~2000。
%                     声阻抗分布按最新要求不再输出。）
%
%    3_signal\      OACA_signal.fig      中心单阵元缺陷回波信号（处理前后同一张图：
%                   OACA_signal_data.mat   处理前黑色虚线、处理后红色实线）
%                  original_signal_data.mat
%                   （两张波形数据 mat：<tag>_signal_data.mat）
%                  original_full_data.mat   相控阵数据（与 simdata_*.mat 同格式：
%                  OACA_full_data.mat       rf / sim_config / medium / input_signal /
%                                           t_axis / x_axis / y_axis）
%
%    4_Image\        original_image.fig   original_image_data.mat
%                   OACA_image.fig       OACA_image_data.mat
%                   original_axial_curve.fig  OACA_axial_curve.fig
%                   （曲线数值同时存于对应 *_image_data.mat / *_axial_curve.mat）
%
%    5_video\       original_propagation_movie.avi
%                   OACA_propagation_movie.avi
%
%  目录名集中定义在第 2 节 CFG.folder，如需改名只改那一处。
%
%
% =========================================================================

clc;
clear;
close all;

%% ---------------- 0. 代码目录与工程根目录 ----------------
%  本文件放在 1_code 里，工程根目录是它的上一级（2_skull_model 等输出文件夹所在处）。
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
code_dir = script_dir;                                   % 代码目录（1_code）
root_dir = fileparts(script_dir);                        % 工程根目录
if isempty(root_dir) || ~exist(root_dir, 'dir')
    root_dir = script_dir;                               % 兜底：代码就放在根目录里
end

%% =====================================================================
%% 1. precheck 按钮
%% =====================================================================
%  0 = 仅查看（只导出模型声速/密度分布图，不跑仿真）
%  1 = 跑完全仿真（仿真 + 相控阵数据 + 动画 + 后处理）
CFG.precheck = 1;

%  是否保存声场传播动画 avi（仅 precheck = 1 时生效）
%  1 = 保存，0 = 不保存。动画文件较大、录制较慢，不需要时置 0。
CFG.save_animation = 1;

%  需要运行的处理状态（可只保留其中一个，例如 {'original'}）
CFG.run_cases = {'original', 'OACA'};

%% ---------------- 2. 目录与文件命名 ----------------
CFG.code_dir  = code_dir;                        % 代码目录（k-Wave、配置快照在这里）
CFG.root_dir  = root_dir;                        % 工程根目录（4 个输出文件夹的上一级）
CFG.kwave_dir = fullfile(code_dir, 'k-wave-toolbox-version-1.4', 'k-Wave');

CFG.folder = struct( ...
    'model',  '2_skull_model', ...               % 文件夹 1：介质分布模型图
    'signal', '3_signal', ...                    % 文件夹 2：单阵元信号 + 相控阵数据
    'image',  '4_Image', ...                     % 文件夹 3：成像图与曲线
    'video',  '5_video');                        % 文件夹 4：声场传播动画
%  说明：目录名沿用现有工程目录，要改名只改上面 image 一项即可，
%        其它地方都从这里取，不会漏改。

%% ---------------- 3. 介质声学参数（OACA 处理前 / 处理后）----------------
%  OACA 处理前：颅骨声速 3110 m/s，密度 1975 kg/m^3
%  OACA 处理后：颅骨声速 1680 m/s，密度 1152 kg/m^3
%  其余两层（耦合凝胶、大脑）在两套模型中保持不变。
%  注意：下面 struct 里的数值就是唯一数据来源，改这里即可，注释同步改。
CFG.media = struct( ...
    'c_gel',   1490,    'rho_gel',   998, ...     % 耦合凝胶  c[m/s], rho[kg/m^3]
    'c_brain', 1540,    'rho_brain', 1050);       % 大脑      c[m/s], rho[kg/m^3]

CFG.cases = struct( ...
    'tag',       {'original',            'OACA'}, ...
    'label',     {'OACA 处理前',         'OACA 处理后'}, ...
    'short',     {'original',            'OACA'}, ...
    'en',        {'before',              'after'}, ...  % 图标题用（Speed before/after OACA）
    'c_skull',   {3110,                1680}, ...
    'rho_skull', {1975,                1152});

%% ---------------- 4. 介质几何参数 ----------------
%  自上而下的分层（曲面以阵列中心为原点，R 为曲率半径）：
%    凝胶   ：上表面平面，下表面球面 R = 5.00 mm
%    颅骨   ：上表面 R = 5.00 mm，下表面 R = 4.75 mm（厚度 0.25 mm）
%    大脑   ：上表面 R = 4.75 mm，下表面为计算域平面底边
CFG.geometry = struct( ...
    'gel_thickness',      1.00e-3, ...   % 凝胶在中心处的厚度 [m]
    'skull_thickness',    0.25e-3, ...   % 颅骨厚度 [m]
    'brain_thickness',    5.00e-3, ...   % 大脑在中心处的厚度 [m]
    'radius_gel_skull',   5.00e-3, ...   % 凝胶/颅骨界面曲率半径 [m]
    'radius_skull_brain', 4.75e-3);      % 颅骨/大脑界面曲率半径 [m]

%% ---------------- 5. 相控阵配置 ----------------
CFG.array = struct( ...
    'n_elements',      32, ...           % 阵元数
    'pitch',           50e-6, ...        % 阵元中心间距 [m]
    'center_freq',     15e6, ...         % 中心频率 [Hz]
    'n_cycles',        1, ...            % 激励周期数（单周期脉冲）
    'source_pressure', 1.0e3, ...        % 单阵元激励峰值声压 [Pa]
    'element_width',   50e-6);           % 阵元宽度（点源近似时仅用于说明）[m]
%  说明：仿真中每个阵元用 1 个网格点作为点源 / 点接收，
%        因此阵元位置直接由 pitch 决定，element_width 不参与计算。

%% ---------------- 6. 缺陷（散射体）分布 ----------------
%  每个缺陷是一个球形散射体（二维仿真中为圆截面），材料参数独立设置。
%  深度 depth 以阵列表面为起点，横向 lateral 相对阵列中心。
%  需要多个缺陷时，按下面注释的方式追加即可（也可用于模拟不同缺陷分布）。
CFG.defects = struct( ...
    'depth',    5.0e-3, ...              % 缺陷中心深度 [m]（阵列正下方 5 mm）
    'lateral',  0.0, ...                 % 缺陷中心横向位置 [m]（0 = 阵列中心）
    'diameter', 0.1e-3, ...              % 缺陷直径 [m]
    'c',        2000, ...                % 缺陷声速 [m/s]
    'rho',      1500, ...                % 缺陷密度 [kg/m^3]
    'label',    '散射体');
%
%  追加第 2 个缺陷的写法（取消注释即可）：
%  CFG.defects(2) = struct('depth', 3.0e-3, 'lateral', 0.6e-3, ...
%      'diameter', 0.1e-3, 'c', 2000, 'rho', 1500, 'label', '散射体 2');
%
%  缺陷分布检查结果（与颅骨/脑层是否重叠）会在运行时逐条打印。

%% ---------------- 7. 计算网格与时间步 ----------------
%  dt（k-Wave 时间步进）仍按 CFL 自动算：dt = dx / min(c) / dt_div ≈ 0.2237 ns
%  （即内部采样率 ≈ 4.47 GHz），保证数值稳定；这个步进不变。
%  导出的相控阵数据在保存前重采样到下面的 fs_out（300 MHz），
%  所以数据文件里的 rf / t_axis 是 300 MHz、0~8 us 的等间隔采样。
CFG.grid = struct( ...
    'dx',             1.0e-6, ...        % 网格尺寸 [m]（深度方向）
    'lateral_width',  3.002e-3, ...      % 横向计算宽度 [m]（3002 个网格点）
    'dt_div',         3, ...             % dt = dx / min(c) / dt_div（仿真步进，不变）
    't_end',          8e-6, ...          % 记录时长 [s]（采样时间 8 us）
    'fs_out',         300e6, ...         % 导出数据的采样频率 [Hz]
    'pml_size',       50, ...            % PML 厚度 [网格点]
    'pml_alpha',      1e3);              % PML 吸收系数
%  说明：fs_out 只影响保存的 rf 与后处理；k-Wave 内部仍用 dt 逐点推进。
%        8 us 对应单程最深 c*t/2 = 1540*8e-6/2 = 6.16 mm，成像 z_range 取 6.0 mm 合适。

%% ---------------- 8. 仿真选项 ----------------
CFG.sim = struct( ...
    'use_gpu',          true, ...        % 有可用 GPU 时自动使用
    'save_checkpoint',  true, ...        % 逐阵元断点续算（长仿真强烈建议开）
    'checkpoint_every', 8, ...           % 每 N 个阵元写一次断点
    'resume',           true, ...        % 存在断点时从断点继续
    'save_full_data',   true, ...        % 保存相控阵数据（3_signal\*_full_data.mat）
    'verbose',          true);           % 打印进度

%% ---------------- 9. 声场传播动画 ----------------
%  动画由一次额外的仿真录制：默认只激励中心阵元（相当于一次单阵元发射的点源
%  扩散），也可以改成 'full' 让 32 个阵元同时同相激励（全孔径平面波）。
%  两种方式都能看到波前穿过凝胶 -> 颅骨 -> 大脑的过程，以及颅骨像差的影响。
%
%  坐标格式：本地 k-Wave 副本已改过动画绘图（kspaceFirstOrder2D.m 里带
%  "[本工程修改]" 标记），横轴 x / mm、纵轴 z / mm 且深度从阵列表面 0 起，
%  字体 Times New Roman 12，与 2_skull_model 的模型图一致。
%  若日后整体替换 k-Wave 工具箱，这段改动会丢失，动画会退回默认的
%  "x-position / y-position" 坐标，按标记再把那几行改回来即可。
CFG.animation = struct( ...
    'aperture',   'center', ...            % 'full' = 全孔径同时激励；'center' = 只激励中心阵元
    'plot_freq',  60, ...                % 每 N 个时间步录一帧（越小越平滑、文件越大）
    'frame_rate', 25, ...                % 动画帧率 [fps]
    'plot_scale', 'auto', ...            % 显示色标：'auto' 或 [min max]（Pa）
    'profile',    'Motion JPEG AVI', ... % VideoWriter 配置
    'figure_position', [100 50 560 860]);% 动画图窗位置与大小（竖长视场）

%% ---------------- 10. 后处理 / 成像参数 ----------------
%  相控阵检测参数与第 5 节共用；这里只补充成像专用设置。
CFG.imaging = struct( ...
    'x_range',     [-1.5e-3, 1.5e-3], ...% 成像横向范围 [m]
    'z_range',     [0, 6.0e-3], ...      % 成像深度范围 [m]
    'step',        0.02e-3, ...          % 成像网格步长 [m]（手动设置）
    'db_range',    40, ...               % 图像显示动态范围 [dB]
    'colormap',    'jet', ...            % 图像色图
    'center_element', [], ...            % 中心单阵元序号（[] = 自动取第 n/2 个）
    'axial_lateral', 0.0);               % 纵向幅值分布所在的横向位置 [m]
%  说明：
%    * 成像介质声速自动取 CFG.media.c_brain（大脑声速），与检测对象一致；
%    * 采样率自动取仿真数据的时间步进 fs = 1/dt；
%    * 成像网格由 x_range / z_range / step 手动给定；
%    * 深度范围建议不超过计算域（凝胶 1.0 + 颅骨 0.25 + 大脑 5.0 = 6.25 mm），
%      再往下是 PML 吸收层，延时叠加在那里会出现边界伪影，不是真实回波。

%% ---------------- 11. 输出选项 ----------------
CFG.output = struct( ...
    'model_figures',   true, ...         % 导出模型分布图（2_skull_model）
    'export_png',      true, ...         % 除 .fig 外再导出一张同名 .png 预览
    'close_figures',   true);            % 出图后关闭图窗，避免堆积

% ---- 模型分布图版式（与 2_skull_model\OACA_Speed.fig 模板一致）----
CFG.model_plot = struct( ...
    'clim_speed',      [1400 3400], ...  % 声速色标范围（处理前后统一）
    'clim_density',    [900 2000], ...   % 密度色标范围（处理前后统一）
    'clim_impedance',  [1.4e6 6.2e6], ...% 声阻抗色标（仅在 draw_impedance = true 时用）
    'font_name',       'Times New Roman', ...
    'font_size',       12, ...           % 坐标区字号（标题/轴标签自动 1.1 倍）
    'figure_size',     [277.6 404], ...  % 图窗宽高 [px]
    'draw_impedance',  false);           % 声阻抗分布不再出图

%% =====================================================================
%% 12. 检查配置 + 打印工程信息
%% =====================================================================
if ~ismember(CFG.precheck, [0, 1])
    error('precheck 只能取 0（仅查看）或 1（跑完全仿真），当前为 %g。', CFG.precheck);
end

fprintf('\n================ 颅骨超声相控阵仿真工程 ================\n');
fprintf('工程根目录：%s\n', CFG.root_dir);
fprintf('precheck   ：%d（%s）\n', CFG.precheck, ...
    iif_str(CFG.precheck == 1, '跑完全仿真', '仅查看，不跑仿真'));
fprintf('动画保存   ：%d（%s）\n', CFG.save_animation, ...
    iif_str(CFG.save_animation == 1, '保存 avi', '不保存'));
fprintf('运行状态   ：%s\n', strjoin(CFG.run_cases, ' / '));
fprintf('--------------------------------------------------------\n');
fprintf('相控阵：%d 阵元，间距 %.1f um，中心频率 %.1f MHz，%.0f 周期，源压 %.0f Pa\n', ...
    CFG.array.n_elements, CFG.array.pitch * 1e6, CFG.array.center_freq / 1e6, ...
    CFG.array.n_cycles, CFG.array.source_pressure);
fprintf('孔径   ：%.2f mm（32 阵元）\n', ...
    (CFG.array.n_elements - 1) * CFG.array.pitch * 1e3);
for ic = 1:numel(CFG.cases)
    fprintf('介质 %-9s：凝胶 %.0f/%.0f，颅骨 %.1f/%.0f，大脑 %.0f/%.0f  (c/rho)\n', ...
        CFG.cases(ic).short, CFG.media.c_gel, CFG.media.rho_gel, ...
        CFG.cases(ic).c_skull, CFG.cases(ic).rho_skull, ...
        CFG.media.c_brain, CFG.media.rho_brain);
end
fprintf('缺陷   ：%d 个\n', numel(CFG.defects));
for id = 1:numel(CFG.defects)
    d = CFG.defects(id);
    fprintf('   [%d] %s：深度 %.2f mm，横向 %.2f mm，直径 %.3f mm，c = %.0f m/s，rho = %.0f kg/m^3\n', ...
        id, d.label, d.depth * 1e3, d.lateral * 1e3, d.diameter * 1e3, d.c, d.rho);
end
fprintf('成像   ：x %s mm，z %s mm，步长 %.3f mm，介质声速 = 大脑 %.0f m/s\n', ...
    mat2str(CFG.imaging.x_range * 1e3), mat2str(CFG.imaging.z_range * 1e3), ...
    CFG.imaging.step * 1e3, CFG.media.c_brain);
fprintf('========================================================\n\n');

% 建立输出目录（4 个文件夹）
create_output_folders(CFG);

% 保存配置快照：放在代码目录，skull_signal_postprocess.m 单独运行时会自动读取它
cfg_snapshot = fullfile(CFG.code_dir, 'skull_model_cfg.mat');
save(cfg_snapshot, 'CFG');
fprintf('配置快照已保存：%s\n\n', cfg_snapshot);

%% =====================================================================
%% 13. 运行
%% =====================================================================
if CFG.precheck == 1

    % ---- 13.1 完整仿真：建模 + 出图 + k-Wave + 相控阵数据 + 动画 ----
    result = skull_simulation_function(CFG, struct('mode', 'simulate'));

    % ---- 13.2 后处理：单阵元信号 + DAS 成像 + 中心纵向幅值分布 ----
    post = skull_signal_postprocess(CFG);

    fprintf('\n================ 仿真 + 后处理全部完成 ================\n');
    fprintf('模型分布图：%s\n', result.folder.model);
    fprintf('相控阵数据：%s\n', result.folder.signal);
    fprintf('声场动画  ：%s\n', result.folder.video);
    fprintf('信号图    ：%s\n', post.folder.signal);
    fprintf('成像/曲线 ：%s\n', post.folder.image);

else

    % ---- precheck = 0：仅建模 + 出模型分布图，不跑仿真 ----
    result = skull_simulation_function(CFG, struct('mode', 'model'));

    fprintf('\n================ precheck = 0：仅查看，未跑仿真 ================\n');
    fprintf('已导出模型分布图：%s\n', result.folder.model);
    fprintf('如需跑完全仿真，请把本文件第 1 节的 CFG.precheck 改为 1 后重新运行。\n');
    fprintf('如需单独检查后处理效果（用现有数据），可先运行：\n');
    fprintf('    skull_signal_postprocess\n');
    fprintf('它会优先用 3_signal 下的 *_full_data.mat，没有则回退到根目录的 simdata_*.mat。\n');
end

%% =====================================================================
%% 局部函数
%% =====================================================================
function create_output_folders(CFG)
%CREATE_OUTPUT_FOLDERS 建立 4 个输出文件夹
names = {CFG.folder.model, CFG.folder.signal, CFG.folder.image, CFG.folder.video};
for i = 1:numel(names)
    d = fullfile(CFG.root_dir, names{i});
    if ~exist(d, 'dir')
        mkdir(d);
        fprintf('已建立目录：%s\n', d);
    end
end
end


function out = iif_str(cond, a, b)
%IIF_STR 条件为真返回 a，否则返回 b
if cond
    out = a;
else
    out = b;
end
end

