function layer_params_opt = fit_layer_parameters_LM(layer_params_init, exp_data_multi, input_params)
% FIT_LAYER_PARAMETERS_LM 使用LM算法拟合层参数（SiO2分解版本）

    fprintf('=== 开始层参数拟合 (SiO2分解模式) ===\n');
    fprintf('优化参数数量: %d\n', length(layer_params_init));
    fprintf('实验数据组数: %d\n', length(exp_data_multi));
    
    % 设置优化选项
    options = optimoptions('lsqnonlin', ...
                         'Display', 'iter', ...
                         'MaxIterations', 60, ...
                         'MaxFunctionEvaluations', 100000, ...
                         'FunctionTolerance', 1e-12, ...
                         'StepTolerance', 1e-12, ...
                         'Algorithm', 'levenberg-marquardt');
    
    % 设置参数边界
    [lb, ub] = set_layer_parameter_bounds_sio2(layer_params_init, input_params);
    
    % 定义目标函数
    objective = @(params) layer_params_objective(params, exp_data_multi, input_params);
    
    % 执行优化
    [layer_params_opt, resnorm, ~, exitflag, output] = ...
        lsqnonlin(objective, layer_params_init, lb, ub, options);
    
    % 输出结果
    fprintf('\n=== 优化完成 ===\n');
    fprintf('退出标志: %d\n', exitflag);
    fprintf('函数评估次数: %d\n', output.funcCount);
    fprintf('最终残差平方和: %.6e\n', resnorm);
    
    % 显示参数对比
    display_layer_optimization_results_sio2(layer_params_init, layer_params_opt, input_params);
end

function [lb, ub] = set_layer_parameter_bounds_sio2(layer_params_init, input_params)
% 设置层参数边界（SiO2分解版本）
    lb = zeros(size(layer_params_init));
    ub = inf(size(layer_params_init));
    
    param_idx = 1;
    for layer = 1:input_params.num_layers
        layer_elements = input_params.layer_info{layer}.elements;
        
        % 检查是否包含SiO2
        if ismember('SiO2', layer_elements)
            % === SiO2层边界设置 ===
            
            % SiO2分解密度边界
            sio2_bulk_density = input_params.bulk_density('SiO2');
            si_bulk_density = input_params.bulk_density('Si');
            o_bulk_density = input_params.bulk_density('O');
            
            % Si密度边界
            lb(param_idx) = si_bulk_density * 0.1;
            ub(param_idx) = sio2_bulk_density * 1.5;  % 允许超过纯Si密度
            
            % 共享σL边界
            lb(param_idx+1) = 0.1;
            ub(param_idx+1) = 10;
            
            % Si σR边界
            lb(param_idx+2) = 0.1;
            ub(param_idx+2) = 10;
            
            % 共享z0边界
            z0_init = input_params.initial_z0(layer);
            layer_thickness = input_params.layer_info{layer}.thickness;
            lb(param_idx+3) = z0_init - layer_thickness;
            ub(param_idx+3) = z0_init + layer_thickness;
            
            % 共享厚度边界
            lb(param_idx+4) = max(0.1, layer_thickness * 0.2);
            ub(param_idx+4) = layer_thickness * 1.5;
            
            % O密度边界
            lb(param_idx+5) = o_bulk_density * 0.1;
            ub(param_idx+5) = sio2_bulk_density * 1.5;
            
            % O σR边界
            lb(param_idx+6) = 0.1;
            ub(param_idx+6) = 10;
            
            param_idx = param_idx + 7;
            
            % 处理该层其他元素
            other_elements = setdiff(layer_elements, {'SiO2'});
            for elem_idx = 1:length(other_elements)
                elem = other_elements{elem_idx};
                bulk_dens = input_params.bulk_density(elem);
                z0_init = input_params.initial_z0(layer);
                layer_thickness = input_params.layer_info{layer}.thickness;
                
                % 普通元素4参数边界
                lb(param_idx:param_idx+3) = [bulk_dens*0.1, 0.1, 0.1, z0_init-layer_thickness];
                ub(param_idx:param_idx+3) = [bulk_dens*10, 10, 10, z0_init+layer_thickness];
                
                param_idx = param_idx + 4;
            end
            
        else
            % === 普通层边界设置 ===
            for elem_idx = 1:length(layer_elements)
                elem = layer_elements{elem_idx};
                bulk_dens = input_params.bulk_density(elem);
                z0_init = input_params.initial_z0(layer);
                layer_thickness = input_params.layer_info{layer}.thickness;
                
                % 普通元素4参数边界
                lb(param_idx:param_idx+3) = [bulk_dens*0.1, 0.1, 0.1, z0_init-layer_thickness];
                ub(param_idx:param_idx+3) = [bulk_dens*10, 10, 10, z0_init+layer_thickness];
                
                param_idx = param_idx + 4;
            end
        end
    end
    
    % 基底σL边界
    if length(layer_params_init) >= param_idx
        lb(param_idx) = 0.1;
        ub(param_idx) = 10.0;
    end
end

function display_layer_optimization_results_sio2(params_init, params_opt, input_params)
% 显示层参数优化结果对比（SiO2分解版本）
    fprintf('\n=== 层参数优化结果对比 (SiO2分解模式) ===\n');
    
    param_idx = 1;
    for layer = 1:input_params.num_layers
        layer_elements = input_params.layer_info{layer}.elements;
        
        % 检查是否包含SiO2
        if ismember('SiO2', layer_elements)
            fprintf('\n第%d层 SiO2 (分解为Si+O):\n', layer);
            
            % Si密度
            change_percent = (params_opt(param_idx) - params_init(param_idx)) / params_init(param_idx) * 100;
            fprintf('  Si ρ     : %.4f → %.4f (变化: %+.1f%%)\n', ...
                    params_init(param_idx), params_opt(param_idx), change_percent);
            
            % 共享σL
            change_percent = (params_opt(param_idx+1) - params_init(param_idx+1)) / params_init(param_idx+1) * 100;
            fprintf('  σL(共享) : %.4f → %.4f (变化: %+.1f%%)\n', ...
                    params_init(param_idx+1), params_opt(param_idx+1), change_percent);
            
            % Si σR
            change_percent = (params_opt(param_idx+2) - params_init(param_idx+2)) / params_init(param_idx+2) * 100;
            fprintf('  Si σR    : %.4f → %.4f (变化: %+.1f%%)\n', ...
                    params_init(param_idx+2), params_opt(param_idx+2), change_percent);
            
            % 共享z0
            change_percent = (params_opt(param_idx+3) - params_init(param_idx+3)) / params_init(param_idx+3) * 100;
            fprintf('  z₀(共享) : %.4f → %.4f (变化: %+.1f%%)\n', ...
                    params_init(param_idx+3), params_opt(param_idx+3), change_percent);
            
            % 共享厚度
            change_percent = (params_opt(param_idx+4) - params_init(param_idx+4)) / params_init(param_idx+4) * 100;
            fprintf('  t(共享)  : %.4f → %.4f (变化: %+.1f%%)\n', ...
                    params_init(param_idx+4), params_opt(param_idx+4), change_percent);
            
            % O密度
            change_percent = (params_opt(param_idx+5) - params_init(param_idx+5)) / params_init(param_idx+5) * 100;
            fprintf('  O ρ      : %.4f → %.4f (变化: %+.1f%%)\n', ...
                    params_init(param_idx+5), params_opt(param_idx+5), change_percent);
            
            % O σR
            change_percent = (params_opt(param_idx+6) - params_init(param_idx+6)) / params_init(param_idx+6) * 100;
            fprintf('  O σR     : %.4f → %.4f (变化: %+.1f%%)\n', ...
                    params_init(param_idx+6), params_opt(param_idx+6), change_percent);
            
            param_idx = param_idx + 7;
            
            % 处理该层其他元素
            other_elements = setdiff(layer_elements, {'SiO2'});
            for elem_idx = 1:length(other_elements)
                elem = other_elements{elem_idx};
                fprintf('第%d层 %s:\n', layer, elem);
                
                % 显示4个参数
                param_names = {'ρ', 'σL', 'σR', 'z₀'};
                for p = 1:4
                    change_percent = (params_opt(param_idx) - params_init(param_idx)) / params_init(param_idx) * 100;
                    fprintf('  %-8s: %.4f → %.4f (变化: %+.1f%%)\n', ...
                            param_names{p}, params_init(param_idx), params_opt(param_idx), change_percent);
                    param_idx = param_idx + 1;
                end
            end
            
        else
            % 普通层处理
            for elem_idx = 1:length(layer_elements)
                elem = layer_elements{elem_idx};
                fprintf('\n第%d层 %s:\n', layer, elem);
                
                param_names = {'ρ', 'σL', 'σR', 'z₀'};
                for p = 1:4
                    change_percent = (params_opt(param_idx) - params_init(param_idx)) / params_init(param_idx) * 100;
                    fprintf('  %-8s: %.4f → %.4f (变化: %+.1f%%)\n', ...
                            param_names{p}, params_init(param_idx), params_opt(param_idx), change_percent);
                    param_idx = param_idx + 1;
                end
            end
        end
    end
    
    % 显示基底参数
    if length(params_init) >= param_idx
        fprintf('\n基底参数:\n');
        change_percent = (params_opt(param_idx) - params_init(param_idx)) / params_init(param_idx) * 100;
        fprintf('  σL       : %.4f → %.4f (变化: %+.1f%%)\n', ...
                params_init(param_idx), params_opt(param_idx), change_percent);
    end
end