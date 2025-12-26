function main_layer_parameter_fitting()
% MAIN_LAYER_PARAMETER_FITTING 层参数拟合主程序（SiO2分解+平滑显示）

    clc; clear; close all;
    
    fprintf('===============================================\n');
    fprintf('  层参数X射线反射率拟合程序\n');
    fprintf('  (SiO2分解+平滑显示版本)\n');
    fprintf('===============================================\n');
    
    try
        %% 步骤1：参数输入
        fprintf('\n【步骤1/5】参数输入\n');
        input_params = input_parameters_layer();
        
        %% 步骤2：构造初始参数向量（SiO2分解）
        fprintf('\n【步骤2/5】构造初始参数向量 (SiO2分解模式)\n');
        [layer_params_init, ~] = construct_initial_layer_params(input_params);
        
        %% 步骤3：加载实验数据
        fprintf('\n【步骤3/5】加载多波长实验数据\n');
        exp_data_multi = load_experimental_data_multi(input_params);
        
        %% 步骤4：执行拟合（SiO2分解）
        fprintf('\n【步骤4/5】开始层参数拟合 (SiO2分解模式)\n');
        layer_params_opt = fit_layer_parameters_LM(layer_params_init, exp_data_multi, input_params);
        
        %% 步骤5：生成完整结果报告（SiO2分解+平滑显示）
        fprintf('\n【步骤5/5】生成完整结果报告 (SiO2分解+平滑显示)\n');
        
        % 构建初始和优化后的密度分布
        grid_model_init = build_density_from_layer_params(layer_params_init, input_params);
        grid_model_opt = build_density_from_layer_params(layer_params_opt, input_params);
        
        % 调用SiO2分解+平滑显示版本的输出函数
        generate_fitting_results_with_MF_R(grid_model_init, grid_model_opt, exp_data_multi, input_params);

        fprintf('\n===============================================\n');
        fprintf('  层参数拟合程序运行完成！\n');
        fprintf('  (SiO2分解+平滑显示版本)\n');
        fprintf('===============================================\n');
        
    catch ME
        fprintf('\n❌ 程序执行出错：%s\n', ME.message);
        fprintf('错误位置：%s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
        if length(ME.stack) > 1
            fprintf('调用链：%s (第%d行)\n', ME.stack(2).name, ME.stack(2).line);
        end
        rethrow(ME);
    end
end