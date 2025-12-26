function input_params = input_parameters_layer()
% INPUT_PARAMETERS_LAYER 层参数拟合的参数输入函数
    
    fprintf('=== 层参数拟合 - 参数输入 ===\n');
    
    %% 物理常数
    re = 2.82e-13;  % 经典电子半径 (cm)
    NA = 6.02214076e23;  % 阿伏伽德罗常数 (atoms/mol)
    
    %% 能量和波长对应关系
    energies = [90, 220, 410, 640, 1125, 1625];  % eV
    wavelengths = [13.78, 5.636, 3.024, 1.937, 1.102, 0.763];  % nm
    fprintf('使用%d个波长进行多波长拟合\n', length(wavelengths));
    
    %% ASF数据表（f1, f2对应6个能量点）
    asf_data = struct();
    % Si
    asf_data.Si.f1 = [0.581, 9.63, 13.12, 13.36, 12.78, 11.66];
    asf_data.Si.f2 = [0.436, 8.56, 4.43, 2.26, 0.859, 0.461];
    % Co
    asf_data.Co.f1 = [8.83, 15.74, 17.31, 15.03, 23.71, 26.62];
    asf_data.Co.f2 = [8.89, 7.13, 4.21, 2.61, 11.56, 6.81];
    % Mg
    asf_data.Mg.f1 = [2.24, 10.04, 11.07, 10.83, 9.54, 11.19];
    asf_data.Mg.f2 = [8.11, 5.25, 2.29, 1.11, 0.333, 3.23];
    % O
    asf_data.O.f1 = [5.93, 6.57, 5.81, 6.62, 8.32, 8.36];
    asf_data.O.f2 = [2.81, 0.919, 0.352, 3.49, 1.45, 0.772];
    % C
    asf_data.C.f1 = [4.26, 3.59, 5.59, 6.28, 6.31, 6.22];
    asf_data.C.f2 = [0.784, 0.245, 2.62, 1.34, 0.512, 0.257];
    
    %% 原子质量和块体密度数据库（添加SiO2）
    atomic_mass = containers.Map({'Si', 'Co', 'Mg', 'O', 'C'}, ...
                                {28.085, 58.933, 24.305, 15.999, 12.011});  % g/mol
    bulk_density = containers.Map({'Si', 'Co', 'Mg', 'O', 'C', 'SiO2'}, ...
                                 {2.33, 8.86, 1.74, 1.43, 2.27, 2.65});  % g/cm³
    
    %% 基底参数（固定）
    substrate_density = 2.33;  % 单晶硅密度 g/cm³
    
    %% 网格和扩散参数（先输入，因为需要用于计算）
    fprintf('\n=== 网格和扩散参数 ===\n');
    dz = input('网格步长(nm, 默认0.1): ');
    if isempty(dz), dz = 0.1; end
    
    initial_sigma = input('初始扩散宽度(nm, 默认0.3): ');
    if isempty(initial_sigma), initial_sigma = 0.3; end
    
    buffer = input('缓冲区大小(nm, 默认2.0): ');
    if isempty(buffer), buffer = 2.0; end
    
    % 添加基底σL参数输入
    substrate_sigma_L = input('基底左界面粗糙度(nm, 默认0.3): ');
    if isempty(substrate_sigma_L), substrate_sigma_L = 0.3; end
    
    %% 输入层结构
    fprintf('\n=== 层结构输入 ===\n');
    fprintf('可用薄膜材料: Si, Co, Mg, O, C, SiO2\n');
    fprintf('基底固定为单晶硅(密度=%.2f g/cm³, σL=%.2f nm)\n', substrate_density, substrate_sigma_L);
    fprintf('空气层: 0 - %.1f nm\n', buffer);
    
    num_layers = input('薄膜层数: ');
    layer_info = cell(num_layers, 1);
    initial_z0 = zeros(num_layers, 1);
    
    % 从薄膜起始位置开始计算（不包括buffer）
    cumulative_depth = 0;
    
    for i = 1:num_layers
        fprintf('\n--- 第 %d 层 ---\n', i);
        
        % 输入材料（可以是多种元素）
        while true
            materials_str = input(sprintf('第%d层材料(空格分隔): ', i), 's');
            materials = strsplit(materials_str);
            valid = true;
            for j = 1:length(materials)
                if ~isKey(bulk_density, materials{j})
                    fprintf('材料 "%s" 不可用\n', materials{j});
                    valid = false;
                    break;
                end
            end
            if valid
                layer_info{i}.elements = materials;
                fprintf('选择: %s\n', strjoin(materials, '+'));
                break;
            end
        end
        
        % 输入厚度
        thickness = input(sprintf('第%d层厚度(nm): ', i));
        layer_info{i}.thickness = thickness;
        
        % 计算层中心位置（相对于薄膜起始位置）
        layer_center_relative = cumulative_depth + thickness/2;
        % 加上buffer得到绝对位置
        initial_z0(i) = layer_center_relative + buffer;
        cumulative_depth = cumulative_depth + thickness;
        
        % 显示层的绝对位置范围
        layer_start_abs = cumulative_depth - thickness + buffer;
        layer_end_abs = cumulative_depth + buffer;
        
        fprintf('第%d层设置: %s, %.1f - %.1f nm (厚度%.1f nm, 中心z₀=%.1f nm)\n', ...
                i, strjoin(materials, '+'), layer_start_abs, layer_end_abs, thickness, initial_z0(i));
    end
    
    %% 输出结构
    input_params = struct();
    input_params.re = re;
    input_params.NA = NA;
    input_params.energies = energies;
    input_params.wavelengths = wavelengths;
    input_params.asf_data = asf_data;
    input_params.atomic_mass = atomic_mass;
    input_params.bulk_density = bulk_density;
    input_params.substrate_density = substrate_density;
    input_params.substrate_sigma_L = substrate_sigma_L;  % 添加基底σL
    
    input_params.num_layers = num_layers;
    input_params.layer_info = layer_info;
    input_params.initial_z0 = initial_z0;
    input_params.dz = dz;
    input_params.initial_sigma = initial_sigma;
    input_params.buffer = buffer;
    
    %% 最终确认显示
    fprintf('\n=== 参数确认 ===\n');
    fprintf('空气层: 0 - %.1f nm (密度=0，参与拟合)\n', buffer);
    fprintf('薄膜层结构:\n');
    
    cumulative_display = buffer;
    for i = 1:num_layers
        layer_start = cumulative_display;
        layer_end = cumulative_display + layer_info{i}.thickness;
        fprintf('  第%d层: %s, %.1f - %.1f nm (厚度%.1f nm, 初始z₀=%.1f nm)\n', ...
                i, strjoin(layer_info{i}.elements, '+'), ...
                layer_start, layer_end, layer_info{i}.thickness, initial_z0(i));
        cumulative_display = layer_end;
    end
    
    fprintf('基底: 单晶硅, σL=%.2f nm (参与拟合)\n', substrate_sigma_L);
    fprintf('薄膜总厚度: %.1f nm\n', cumulative_depth);
    fprintf('总计算深度: %.1f nm (空气层 + 薄膜)\n', buffer + cumulative_depth);
    fprintf('网格步长: %.3f nm\n', dz);
    fprintf('初始扩散宽度: %.2f nm\n', initial_sigma);
    fprintf('使用%d个波长进行拟合\n', length(wavelengths));
    
    fprintf('层参数输入完成！\n');
end