function [layer_params_init, param_info] = construct_initial_layer_params(input_params)
% CONSTRUCT_INITIAL_LAYER_PARAMS 构造初始层参数向量（SiO2自动分解）
% 修改：SiO2分解为Si和O，共享σL/z0/thickness，独立ρ/σR

    fprintf('=== 构造初始层参数向量 (SiO2分解模式) ===\n');
    
    layer_params_init = [];
    param_info = {};
    param_count = 0;
    
    for layer = 1:input_params.num_layers
        layer_elements = input_params.layer_info{layer}.elements;
        layer_thickness = input_params.layer_info{layer}.thickness;
        
        fprintf('第%d层 (%s, %.1f nm):\n', layer, strjoin(layer_elements, '+'), layer_thickness);
        
        % 检查是否包含SiO2
        if ismember('SiO2', layer_elements)
            % === SiO2特殊处理：分解为Si和O ===
            fprintf('  检测到SiO2，分解为Si和O (共享σL/z0/thickness)\n');
            
            % SiO2分解密度计算
            sio2_bulk_density = input_params.bulk_density('SiO2'); % 2.65 g/cm³
            M_Si = input_params.atomic_mass('Si');  % 28.085
            M_O = input_params.atomic_mass('O');    % 15.999
            M_SiO2 = M_Si + 2*M_O;                  % 60.083
            
            % Si和O的质量分数
            w_Si = M_Si / M_SiO2;      % ≈ 0.467
            w_O = 2*M_O / M_SiO2;      % ≈ 0.533
            
            % Si和O的体密度
            rho_Si_from_SiO2 = sio2_bulk_density * w_Si;  % ≈ 1.238 g/cm³
            rho_O_from_SiO2 = sio2_bulk_density * w_O;    % ≈ 1.412 g/cm³
            
            % 共享参数
            sigma_L_shared = input_params.initial_sigma;    % 上界面共享
            z0_shared = input_params.initial_z0(layer);     % 中心位置共享
            thickness_shared = layer_thickness;             % 厚度共享
            
            % Si独立参数
            sigma_R_Si = input_params.initial_sigma;        % Si下界面
            
            % O独立参数  
            sigma_R_O = input_params.initial_sigma;         % O下界面
            
            % === 添加分解后的Si参数 ===
            layer_params_init = [layer_params_init; rho_Si_from_SiO2];  % Si密度
            param_count = param_count + 1;
            param_info{param_count} = sprintf('L%d_Si_rho', layer);
            
            layer_params_init = [layer_params_init; sigma_L_shared];     % 共享σL
            param_count = param_count + 1;
            param_info{param_count} = sprintf('L%d_SiO2_sigmaL', layer);
            
            layer_params_init = [layer_params_init; sigma_R_Si];         % Si独立σR
            param_count = param_count + 1;
            param_info{param_count} = sprintf('L%d_Si_sigmaR', layer);
            
            layer_params_init = [layer_params_init; z0_shared];          % 共享z0
            param_count = param_count + 1;
            param_info{param_count} = sprintf('L%d_SiO2_z0', layer);
            
            layer_params_init = [layer_params_init; thickness_shared];   % 共享thickness
            param_count = param_count + 1;
            param_info{param_count} = sprintf('L%d_SiO2_thickness', layer);
            
            % === 添加分解后的O参数 ===
            layer_params_init = [layer_params_init; rho_O_from_SiO2];   % O密度
            param_count = param_count + 1;
            param_info{param_count} = sprintf('L%d_O_rho', layer);
            
            layer_params_init = [layer_params_init; sigma_R_O];          % O独立σR
            param_count = param_count + 1;
            param_info{param_count} = sprintf('L%d_O_sigmaR', layer);
            
            fprintf('  Si: ρ=%.3f, σL=%.3f(共享), σR=%.3f, z₀=%.3f(共享), t=%.3f(共享)\n', ...
                    rho_Si_from_SiO2, sigma_L_shared, sigma_R_Si, z0_shared, thickness_shared);
            fprintf('  O:  ρ=%.3f, σL=%.3f(共享), σR=%.3f, z₀=%.3f(共享), t=%.3f(共享)\n', ...
                    rho_O_from_SiO2, sigma_L_shared, sigma_R_O, z0_shared, thickness_shared);
            
            % 处理该层的其他元素（如果有）
            other_elements = setdiff(layer_elements, {'SiO2'});
            for elem_idx = 1:length(other_elements)
                elem = other_elements{elem_idx};
                
                rho_init = input_params.bulk_density(elem);
                z0_elem = input_params.initial_z0(layer);
                
                % 普通元素的4个参数
                layer_params_init = [layer_params_init; rho_init];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_rho', layer, elem);
                
                layer_params_init = [layer_params_init; input_params.initial_sigma];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_sigmaL', layer, elem);
                
                layer_params_init = [layer_params_init; input_params.initial_sigma];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_sigmaR', layer, elem);
                
                layer_params_init = [layer_params_init; z0_elem];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_z0', layer, elem);
                
                fprintf('  %s: ρ=%.3f, σL=%.3f, σR=%.3f, z₀=%.3f\n', ...
                        elem, rho_init, input_params.initial_sigma, ...
                        input_params.initial_sigma, z0_elem);
            end
            
        else
            % === 普通层处理 ===
            for elem_idx = 1:length(layer_elements)
                elem = layer_elements{elem_idx};
                
                if isKey(input_params.bulk_density, elem)
                    rho_init = input_params.bulk_density(elem);
                else
                    error('元素 "%s" 不在bulk_density中', elem);
                end
                
                % 普通元素的4个参数
                layer_params_init = [layer_params_init; rho_init];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_rho', layer, elem);
                
                layer_params_init = [layer_params_init; input_params.initial_sigma];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_sigmaL', layer, elem);
                
                layer_params_init = [layer_params_init; input_params.initial_sigma];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_sigmaR', layer, elem);
                
                layer_params_init = [layer_params_init; input_params.initial_z0(layer)];
                param_count = param_count + 1;
                param_info{param_count} = sprintf('L%d_%s_z0', layer, elem);
                
                fprintf('  %s: ρ=%.3f, σL=%.3f, σR=%.3f, z₀=%.3f\n', ...
                        elem, rho_init, input_params.initial_sigma, ...
                        input_params.initial_sigma, input_params.initial_z0(layer));
            end
        end
    end
    
    % 添加基底σL参数
    substrate_sigma_L_init = input_params.substrate_sigma_L;
    layer_params_init = [layer_params_init; substrate_sigma_L_init];
    param_count = param_count + 1;
    param_info{param_count} = 'Substrate_sigmaL';
    
    fprintf('\n基底参数:\n');
    fprintf('  基底σL: %.3f\n', substrate_sigma_L_init);
    
    fprintf('初始参数向量构造完成！共%d个参数\n', param_count);
    fprintf('SiO2分解说明: 每层SiO2 → 7参数(ρ_Si, σL_共享, σR_Si, z₀_共享, t_共享, ρ_O, σR_O)\n');
    input_params.param_info = param_info;
    
    % 保存SiO2分解标记
    input_params.sio2_decomposed = true;
end