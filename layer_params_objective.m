function residuals = layer_params_objective(layer_params, exp_data_multi, input_params)
% LAYER_PARAMS_OBJECTIVE 层参数拟合的目标函数

    try
        % 从层参数构建密度分布
        grid_model = build_density_from_layer_params(layer_params, input_params);
        
        % 硬编码权重（可手动调节）
        weights = [7, 2.6, 2.5, 1, 1, 2];  % 6个波长的权重
        
        % 确保权重数量匹配
        if length(weights) ~= length(exp_data_multi)
            error('权重数组长度与数据组数不匹配');
        end
        
        % 计算所有波长的残差
        all_residuals = [];
        for i = 1:length(exp_data_multi)
            theta_exp = exp_data_multi(i).theta;
            R_exp = exp_data_multi(i).R;
            energy_index = exp_data_multi(i).energy_index;
            weight = weights(i);
            
            % 计算理论反射率（复用参考代码函数）
            R_calc = calculate_reflectivity(grid_model, input_params, energy_index, theta_exp);
            
            % 计算加权残差
            R_residuals = weight * (R_exp - R_calc) ./ R_exp;
            all_residuals = [all_residuals; R_residuals(:)];
        end
        residuals = all_residuals;
        
    catch ME
        fprintf('目标函数计算出错: %s\n', ME.message);
        total_data_points = sum(arrayfun(@(x) length(x.theta), exp_data_multi));
        residuals = 1e6 * ones(total_data_points, 1);
    end
end