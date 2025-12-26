function [chi_real, chi_imag] = density_to_chi(grid_model, input_params, energy_index)
% DENSITY_TO_CHI 将密度分布转换为极化率分布（SiO2分解版本）
% 修改：移除SiO2特殊处理，因为已分解为Si和O基本元素

    %% 提取参数
    density_grid = grid_model.density_grid;
    all_elements = grid_model.all_elements;  % 现在只包含基本元素：Si, Co, Mg, O, C
    num_points = grid_model.num_points;
    
    % 物理常数
    re = input_params.re;           % 经典电子半径 (cm)
    NA = input_params.NA;           % 阿伏伽德罗常数
    
    % 当前能量的参数
    wavelength = input_params.wavelengths(energy_index);  % nm
    lambda_cm = wavelength * 1e-7;       % 转换为cm (标准单位)
    
    asf_data = input_params.asf_data;
    atomic_mass = input_params.atomic_mass;
    
    %% 初始化输出数组
    chi_real = zeros(num_points, 1);
    chi_imag = zeros(num_points, 1);
    
    %% 遍历每个网格点
    for i = 1:num_points
        chi_point_real = 0;
        chi_point_imag = 0;
        
        % 遍历所有基本元素
        for j = 1:length(all_elements)
            elem = all_elements{j};
            
            % 获取该元素在该点的密度 (g/cm³)
            rho_j = density_grid.(elem)(i);
            
            if rho_j > 1e-10  % 只计算有意义的密度
                
                % 获取原子质量 μⱼ (g/mol)
                if isKey(atomic_mass, elem)
                    mu_j = atomic_mass(elem);
                else
                    error('未找到元素 %s 的原子质量', elem);
                end
                
                % 计算数密度 nⱼ (atoms/cm³)
                n_j = rho_j * NA / mu_j;
                
                % 获取该元素在当前能量的ASF
                if isfield(asf_data, elem)
                    f1_j = asf_data.(elem).f1(energy_index);
                    f2_j = asf_data.(elem).f2(energy_index);
                else
                    error('未找到元素 %s 的ASF数据', elem);
                end
                
                % 计算极化率贡献（标准公式）
                prefactor = (re / pi) * lambda_cm^2;
                
                % 累加该元素对χ的贡献
                contribution_real = prefactor * n_j * f1_j;
                contribution_imag = -prefactor * n_j * f2_j;
                
                chi_point_real = chi_point_real + contribution_real;
                chi_point_imag = chi_point_imag + contribution_imag;
            end
        end
        
        % 保存该点的χ值
        chi_real(i) = chi_point_real;
        chi_imag(i) = chi_point_imag;
    end
    
    %% 数据验证
    if any(isnan(chi_real)) || any(isnan(chi_imag))
        warning('density_to_chi: 检测到NaN值，请检查输入数据');
    end
    
    %% 调试信息
    if any(isinf(chi_real)) || any(isinf(chi_imag))
        warning('density_to_chi: 检测到Inf值，可能存在数值问题');
    end
    
    fprintf('极化率计算完成 - 使用SiO2分解模式（基本元素：%s）\n', strjoin(all_elements, ', '));
end