function generate_fitting_results_with_MF_R(grid_model_init, grid_model_opt, exp_data_multi, input_params)
% GENERATE_FITTING_RESULTS_WITH_MF_R 生成完整的拟合结果报告（SiO2分解+平滑显示）

    %% 创建结果文件夹
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    result_folder = fullfile('fittingresult', timestamp);
    if ~exist('fittingresult', 'dir')
        mkdir('fittingresult');
    end
    mkdir(result_folder);
    fprintf('结果保存到: %s\n', result_folder);
    
    %% 计算MF_R 
    fprintf('计算MF_R评估指标...\n');
    MF_R_individual = zeros(length(exp_data_multi), 1);
    sum_squared_residuals = 0;
    total_points = 0;
    
    for i = 1:length(exp_data_multi)
        theta_exp = exp_data_multi(i).theta;
        R_exp = exp_data_multi(i).R;
        energy_index = exp_data_multi(i).energy_index;
        
        R_calc = calculate_reflectivity(grid_model_opt, input_params, energy_index, theta_exp);
        residuals = (R_exp - R_calc) ./ R_exp;
        MF_R_individual(i) = sum(residuals.^2);
        
        sum_squared_residuals = sum_squared_residuals + sum(residuals.^2);
        total_points = total_points + length(theta_exp);
    end
    MF_R_total = sum_squared_residuals / total_points;
    
    fprintf('总体MF_R: %.6f (SiO2分解模式)\n', MF_R_total);
    
    %% 第一部分：反射率拟合对比图
    fprintf('生成反射率拟合图...\n');
    fig1 = figure('Name', '多波长反射率拟合结果 (SiO2分解)', 'Position', [100, 100, 1200, 800]);
    
    reflectivity_data = struct();
    
    for i = 1:length(exp_data_multi)
        subplot(2, 3, i);
        
        theta_exp = exp_data_multi(i).theta;
        R_exp = exp_data_multi(i).R;
        energy_index = exp_data_multi(i).energy_index;
        wavelength = input_params.wavelengths(energy_index);
        energy = input_params.energies(energy_index);
        
        % 计算连续角度的反射率用于绘图
        theta_max = max(theta_exp) * 1.1;
        theta_calc = linspace(0, theta_max, 300);
        R_init_calc = calculate_reflectivity(grid_model_init, input_params, energy_index, theta_calc);
        R_opt_calc = calculate_reflectivity(grid_model_opt, input_params, energy_index, theta_calc);
        R_calc_exp = calculate_reflectivity(grid_model_opt, input_params, energy_index, theta_exp);
        
        % 绘图
        semilogy(theta_calc, R_init_calc, 'r--', 'LineWidth', 1.5, 'DisplayName', '初始计算');
        hold on;
        semilogy(theta_calc, R_opt_calc, 'b-', 'LineWidth', 2, 'DisplayName', '拟合结果');
        semilogy(theta_exp, R_exp, 'ko', 'MarkerSize', 4, 'DisplayName', '实验数据');
        
        xlabel('角度 (度)');
        ylabel('反射率');
        title(sprintf('%.1f eV (%.3f nm)\nMF_R = %.4f', energy, wavelength, MF_R_individual(i)));
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        xlim([0, theta_max]);
        hold off;
        
        % 存储数据
        field_name = sprintf('wl_%d_nm', round(wavelength*1000));
        reflectivity_data.(field_name).theta_exp = theta_exp;
        reflectivity_data.(field_name).R_exp = R_exp;
        reflectivity_data.(field_name).R_calc = R_calc_exp;
        reflectivity_data.(field_name).theta_smooth = theta_calc;
        reflectivity_data.(field_name).R_smooth = R_opt_calc;
        reflectivity_data.(field_name).MF_R = MF_R_individual(i);
        reflectivity_data.(field_name).energy = energy;
        reflectivity_data.(field_name).wavelength = wavelength;
    end
    
    sgtitle(sprintf('多波长反射率拟合结果 (SiO2分解, MF_R = %.6f)', MF_R_total), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % 保存图片和数据
    savefig(fig1, fullfile(result_folder, 'reflectivity_fitting_sio2_decomposed.fig'));
    
    % 导出反射率数据
    fid = fopen(fullfile(result_folder, 'reflectivity_data_sio2_decomposed.txt'), 'w');
    fprintf(fid, '%% 多波长反射率拟合数据 (SiO2分解模式)\n');
    fprintf(fid, '%% 生成时间: %s\n', datestr(now));
    fprintf(fid, '%% 总体MF_R: %.6f\n', MF_R_total);
    fprintf(fid, '%% 说明: SiO2自动分解为Si和O，共享σL/z0/thickness，独立ρ/σR\n\n');
    
    for i = 1:length(exp_data_multi)
        wavelength = input_params.wavelengths(exp_data_multi(i).energy_index);
        energy = input_params.energies(exp_data_multi(i).energy_index);
        fprintf(fid, '%% 波长: %.3f nm, 能量: %.1f eV, MF_R: %.6f\n', ...
                wavelength, energy, MF_R_individual(i));
        fprintf(fid, 'theta_exp\tR_exp\tR_calc\n');
        
        field_name = sprintf('wl_%d_nm', round(wavelength*1000));
        data = reflectivity_data.(field_name);
        for j = 1:length(data.theta_exp)
            fprintf(fid, '%.6f\t%.6e\t%.6e\n', ...
                    data.theta_exp(j), data.R_exp(j), data.R_calc(j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    close(fig1);
    
    %% 第二部分：密度分布对比图（平滑显示）
    fprintf('生成密度分布对比图 (平滑显示)...\n');
    
    z_grid = grid_model_opt.z_grid;
    density_grid_init = grid_model_init.density_grid;
    density_grid_opt = grid_model_opt.density_grid;
    all_elements = {'Si', 'Co', 'Mg', 'O', 'C'};  % 不包含SiO2
    
    % 检测实际存在的元素
    existing_elements = {};
    for i = 1:length(all_elements)
        elem = all_elements{i};
        if isfield(density_grid_opt, elem) && max(density_grid_opt.(elem)) > 1e-6
            existing_elements{end+1} = elem;
        end
    end
    
    if isempty(existing_elements)
        fprintf('警告：未检测到任何有效的元素密度分布\n');
        return;
    end
    
    % 创建密度分布对比图（平滑显示）
    fig2 = figure('Name', '密度分布对比 (SiO2分解, 平滑显示)', 'Position', [200, 200, 800, 600]);
    
    % 颜色映射
    element_colors = containers.Map({'Si', 'Co', 'Mg', 'O', 'C'}, ...
                                   {[1,0,0], [0,0,1], [0,1,0], [1,0.5,0], [0.5,0,1]});
    
    for i = 1:length(existing_elements)
        elem = existing_elements{i};
        density_init = density_grid_init.(elem);
        density_opt = density_grid_opt.(elem);
        
        if isKey(element_colors, elem)
            color = element_colors(elem);
        else
            color = rand(1,3);
        end
        
        % 使用plot显示平滑曲线（不是stairs台阶）
        plot(z_grid, density_init, '--', 'Color', color, 'LineWidth', 1.5, ...
             'DisplayName', sprintf('%s (初始)', elem));
        hold on;
        plot(z_grid, density_opt, '-', 'Color', color, 'LineWidth', 2.5, ...
             'DisplayName', sprintf('%s (拟合)', elem));
    end
    
    xlabel('深度 (nm)');
    ylabel('密度 (g/cm³)');
    title('密度分布对比 (SiO2分解, 平滑显示)');
    legend('Location', 'best');
    grid on;
    xlim([0, max(z_grid)]);
    ylim([0, max(ylim)*1.05]);
    hold off;
    
    % 保存图片
    savefig(fig2, fullfile(result_folder, 'density_distribution_sio2_smooth.fig'));
    
    % 导出密度数据
    fid = fopen(fullfile(result_folder, 'density_distribution_sio2_smooth.txt'), 'w');
    fprintf(fid, '%% 密度分布对比数据 (SiO2分解, 平滑显示)\n');
    fprintf(fid, '%% 生成时间: %s\n', datestr(now));
    fprintf(fid, '%% 说明: SiO2已分解为独立的Si和O密度分布，采用平滑曲线显示\n');
    fprintf(fid, '%% 网格步长: %.3f nm\n', grid_model_opt.dz);
    
    % 写入表头
    fprintf(fid, 'depth(nm)');
    for i = 1:length(existing_elements)
        elem = existing_elements{i};
        fprintf(fid, '\t%s_init(g/cm3)\t%s_fit(g/cm3)', elem, elem);
    end
    fprintf(fid, '\n');
    
    % 写入数据
    for i = 1:length(z_grid)
        fprintf(fid, '%.6f', z_grid(i));
        for j = 1:length(existing_elements)
            elem = existing_elements{j};
            fprintf(fid, '\t%.6f\t%.6f', density_grid_init.(elem)(i), density_grid_opt.(elem)(i));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    close(fig2);
    
    %% 保存标准density_distribution.mat文件
    fprintf('保存标准density_distribution.mat文件 (SiO2分解)...\n');
    
    % 构建标准兼容的数据结构
    density_data = struct();
    density_data.z_grid = z_grid;
    density_data.timestamp = timestamp;
    density_data.grid_step = input_params.dz;
    density_data.calculation_depth = max(z_grid);
    density_data.source = 'layer_parameter_fitting_sio2_decomposed';
    density_data.sio2_decomposed = true;
    
    % 使用标准字段名格式
    all_elements_standard = {'Si', 'Co', 'Mg', 'O', 'C'};
    fprintf('\n元素密度统计 (SiO2分解):\n');
    
    for i = 1:length(all_elements_standard)
        elem = all_elements_standard{i};
        field_name = sprintf('density_%s', elem);
        
        if isfield(density_grid_opt, elem)
            density_data.(field_name) = density_grid_opt.(elem);
            
            % 统计信息
            max_density = max(density_grid_opt.(elem));
            nonzero_points = sum(density_grid_opt.(elem) > 1e-6);
            fprintf('  %s: 最大密度=%.4f g/cm³, 非零点=%d\n', elem, max_density, nonzero_points);
        else
            density_data.(field_name) = zeros(size(z_grid));
            fprintf('  %s: 无分布\n', elem);
        end
    end
    
    % 保存MAT文件
    save(fullfile(result_folder, 'density_distribution.mat'), 'density_data', '-v7.3');
    fprintf('✅ density_distribution.mat 已保存 (SiO2分解模式)\n');
    
    %% 第三部分：极化率/λ²分布图（平滑显示）
    fprintf('生成极化率/λ²分布图 (平滑显示)...\n');
    fig3 = figure('Name', '极化率/λ²分布 (SiO2分解, 平滑)', 'Position', [300, 300, 1000, 600]);
    
    colors = lines(length(input_params.wavelengths));
    chi_data = struct();
    
    for i = 1:length(input_params.wavelengths)
        wavelength = input_params.wavelengths(i);
        energy = input_params.energies(i);
        
        % 计算极化率
        [chi_real, chi_imag] = density_to_chi(grid_model_opt, input_params, i);
        
        % 计算χ/λ²
        chi_real_normalized = chi_real / (wavelength^2);
        chi_imag_normalized = chi_imag / (wavelength^2);
        
        % 存储数据
        field_name = sprintf('wl_%d_nm', round(wavelength*1000));
        chi_data.(field_name).chi_real_norm = chi_real_normalized;
        chi_data.(field_name).chi_imag_norm = chi_imag_normalized;
        chi_data.(field_name).energy = energy;
        chi_data.(field_name).wavelength = wavelength;
        
        % 绘制实部/λ²（平滑显示）
        subplot(1, 2, 1);
        plot(z_grid, chi_real_normalized, 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('%.1f eV', energy));
        hold on;
        
        % 绘制虚部/λ²（取负号，平滑显示）
        subplot(1, 2, 2);
        plot(z_grid, -chi_imag_normalized, 'Color', colors(i,:), 'LineWidth', 2, ...
             'DisplayName', sprintf('%.1f eV', energy));
        hold on;
    end
    
    subplot(1, 2, 1);
    xlabel('深度 (nm)');
    ylabel('Re(χ/λ²) (nm⁻²)');
    title('极化率实部/λ² (SiO2分解, 平滑)');
    legend('Location', 'best');
    grid on;
    
    subplot(1, 2, 2);
    xlabel('深度 (nm)');
    ylabel('-Im(χ/λ²) (nm⁻²)');
    title('极化率虚部/λ² (SiO2分解, 平滑)');
    legend('Location', 'best');
    grid on;
    
    sgtitle('极化率/λ² 分布 (SiO2分解, 平滑显示)', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 保存图片和数据
    savefig(fig3, fullfile(result_folder, 'chi_over_lambda2_sio2_smooth.fig'));
    
    % 导出极化率数据
    fid = fopen(fullfile(result_folder, 'chi_over_lambda2_sio2_smooth.txt'), 'w');
    fprintf(fid, '%% 极化率/λ²分布数据 (SiO2分解, 平滑显示)\n');
    fprintf(fid, '%% 生成时间: %s\n', datestr(now));
    fprintf(fid, '%% 说明: SiO2分解为Si+O后计算的极化率分布，平滑曲线显示\n\n');
    
    % 写入表头
    fprintf(fid, 'depth(nm)');
    for i = 1:length(input_params.wavelengths)
        energy = input_params.energies(i);
        fprintf(fid, '\tRe_chi/λ²_%.1feV\t-Im_chi/λ²_%.1feV', energy, energy);
    end
    fprintf(fid, '\n');
    
    % 写入数据
    for j = 1:length(z_grid)
        fprintf(fid, '%.6f', z_grid(j));
        for i = 1:length(input_params.wavelengths)
            wavelength = input_params.wavelengths(i);
            field_name = sprintf('wl_%d_nm', round(wavelength*1000));
            fprintf(fid, '\t%.6e\t%.6e', ...
                    chi_data.(field_name).chi_real_norm(j), ...
                    -chi_data.(field_name).chi_imag_norm(j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    close(fig3);
    
    %% 生成总结报告
    fid = fopen(fullfile(result_folder, 'fitting_summary_sio2_decomposed.txt'), 'w');
    fprintf(fid, '层参数X射线反射率拟合结果总结 (SiO2分解+平滑显示)\n');
    fprintf(fid, '=========================================================\n');
    fprintf(fid, '生成时间: %s\n', datestr(now));
    fprintf(fid, '结果文件夹: %s\n\n', result_folder);
    
    fprintf(fid, '拟合质量评估 (MF_R):\n');
    fprintf(fid, '计算公式: MF_R = sum((R_exp - R_calc)/R_exp)^2\n');
    fprintf(fid, '总体MF_R: %.6f\n\n', MF_R_total);
    
    fprintf(fid, '各波长MF_R详情:\n');
    for i = 1:length(exp_data_multi)
        wavelength = input_params.wavelengths(exp_data_multi(i).energy_index);
        energy = input_params.energies(exp_data_multi(i).energy_index);
        fprintf(fid, '%.3f nm (%.1f eV): MF_R = %.6f\n', ...
                wavelength, energy, MF_R_individual(i));
    end
    
    fprintf(fid, '\nSiO2分解处理:\n');
    fprintf(fid, '- 输入: SiO2层\n');
    fprintf(fid, '- 分解: Si + O (质量比 %.1f%% : %.1f%%)\n', 46.7, 53.3);
    fprintf(fid, '- 共享参数: σL(上界面), z₀(中心), thickness(厚度)\n');
    fprintf(fid, '- 独立参数: ρ(密度), σR(下界面)\n');
    fprintf(fid, '- 物理意义: 上界面共同，下界面独立过渡\n');
    
    fprintf(fid, '\n检测到的元素: %s\n', strjoin(existing_elements, ', '));
    
    fprintf(fid, '\n显示特性:\n');
    fprintf(fid, '- 密度分布: 平滑曲线显示，SiO2分解为Si+O独立显示\n');
    fprintf(fid, '- 极化率分布: 平滑曲线显示，基于分解后的Si+O计算\n');
    fprintf(fid, '- 反射率拟合: 初始(虚线) vs 拟合(实线)对比\n');
    
    fprintf(fid, '\n生成文件列表:\n');
    fprintf(fid, '- reflectivity_fitting_sio2_decomposed.fig/.txt: 反射率拟合\n');
    fprintf(fid, '- density_distribution_sio2_smooth.fig/.txt: 密度分布(平滑)\n');
    fprintf(fid, '- chi_over_lambda2_sio2_smooth.fig/.txt: 极化率分布(平滑)\n');
    fprintf(fid, '- density_distribution.mat: 标准兼容MAT文件\n');
    fprintf(fid, '- fitting_summary_sio2_decomposed.txt: 本总结报告\n');
    
    fclose(fid);
    
    fprintf('\n=== SiO2分解+平滑显示拟合结果生成完成 ===\n');
    fprintf('总体MF_R: %.6f\n', MF_R_total);
    fprintf('结果已保存到: %s\n', result_folder);
    fprintf('检测到的元素: %s\n', strjoin(existing_elements, ', '));
    fprintf('SiO2处理: 自动分解为Si+O，共享σL/z₀/thickness，独立ρ/σR\n');
    fprintf('显示模式: 平滑曲线（非台阶）\n');
    
    % 显示各波长MF_R
    fprintf('\n各波长MF_R详情:\n');
    for i = 1:length(exp_data_multi)
        wavelength = input_params.wavelengths(exp_data_multi(i).energy_index);
        energy = input_params.energies(exp_data_multi(i).energy_index);
        fprintf('%.3f nm (%.1f eV): MF_R = %.6f\n', ...
                wavelength, energy, MF_R_individual(i));
    end
end