function R = calculate_reflectivity(grid_model, input_params, energy_index, theta_deg)
% CALCULATE_REFLECTIVITY 使用Parratt递推算法计算X射线反射率
%   支持SiO2化合物和小角度数值稳定性

    %% 输入验证
    if nargin < 4
        error('需要4个输入参数：grid_model, input_params, energy_index, theta_deg');
    end
    
    %% 提取参数
    z_grid = grid_model.z_grid;
    num_points = length(z_grid);
    
    % 当前能量参数
    energy = input_params.energies(energy_index);
    wavelength = input_params.wavelengths(energy_index);  % nm
    k0 = 2 * pi / wavelength;  % 波矢 (1/nm)
    
    %% 小角度处理
    theta_deg_safe = theta_deg;
    min_angle = 0.001;
    theta_deg_safe(theta_deg_safe < min_angle) = min_angle;
    
    theta_rad = theta_deg_safe * pi / 180;
    cos_theta = cos(theta_rad);
    
    %% 转换密度分布为极化率
    [chi_real, chi_imag] = density_to_chi(grid_model, input_params, energy_index);
    
    % 检查χ数值
    if any(~isfinite(chi_real)) || any(~isfinite(chi_imag))
        warning('χ存在非有限值，可能导致反射率计算失败');
    end
    
    %% 计算复折射率
    delta = -chi_real / 2;
    beta = -chi_imag / 2;
    n_complex = 1 + delta + 1i * beta;
    
    %% 计算反射率
    num_angles = length(theta_rad);
    R = zeros(num_angles, 1);
    
    for angle_idx = 1:num_angles
        cos_theta_i = cos_theta(angle_idx);
        
        try
            % 计算每层的q_z
            qz = zeros(num_points, 1);
            for j = 1:num_points
                n_j = n_complex(j);
                
                % 数值稳定的q_z计算
                n2_minus_cos2 = n_j^2 - cos_theta_i^2;
                
                if real(n2_minus_cos2) >= 0
                    qz(j) = k0 * sqrt(n2_minus_cos2);
                else
                    qz(j) = k0 * 1i * sqrt(-n2_minus_cos2);
                end
                
                % 确保虚部为正
                if imag(qz(j)) < 0
                    qz(j) = -qz(j);
                end
            end
            
            %% Parratt递推算法
            r = 0;  % 基底反射系数为0
            
            % 从底部向顶部递推
            for j = num_points-1:-1:1
                % 层间界面的Fresnel系数
                qz_j = qz(j);
                qz_j1 = qz(j+1);
                
                if abs(qz_j + qz_j1) < 1e-15
                    r_fresnel = 0;
                else
                    r_fresnel = (qz_j - qz_j1) / (qz_j + qz_j1);
                end
                
                % 层厚度
                if j == 1
                    thickness = z_grid(2) - z_grid(1);
                else
                    thickness = z_grid(j+1) - z_grid(j);
                end
                
                % 相位因子
                phase = exp(2i * qz_j1 * thickness);
                
                % 数值稳定性限制
                if abs(phase) > 1e10
                    phase = sign(phase) * 1e10;
                end
                
                % Parratt递推公式
                denominator = 1 + r * r_fresnel * phase;
                if abs(denominator) < 1e-15
                    r = r_fresnel;
                else
                    r = (r_fresnel + r * phase) / denominator;
                end
                
                if ~isfinite(r)
                    r = r_fresnel;
                    break;
                end
            end
            
            % 计算反射率
            R_current = abs(r)^2;
            
            % 物理约束
            if R_current < 0
                R_current = 0;
            elseif R_current > 1
                R_current = 1;
            end
            
            R(angle_idx) = R_current;
            
        catch
            R(angle_idx) = 0;
        end
    end
    
    %% 小角度插值修正
    if min(theta_deg) < min_angle
        zero_indices = theta_deg < min_angle;
        if any(zero_indices)
            safe_value = R(find(~zero_indices, 1, 'first'));
            R(zero_indices) = safe_value;
        end
    end
    
    %% 最终数值修复
    nan_count = sum(~isfinite(R));
    if nan_count > 0
        valid_indices = isfinite(R);
        if any(valid_indices)
            R(~valid_indices) = interp1(theta_deg(valid_indices), R(valid_indices), ...
                                       theta_deg(~valid_indices), 'linear', 'extrap');
        else
            R(:) = 1e-10;
        end
        warning('calculate_reflectivity: 反射率数值异常，已修复%d个点', nan_count);
    end
end