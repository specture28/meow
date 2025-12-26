function exp_data_multi = load_experimental_data_multi(input_params)
% LOAD_EXPERIMENTAL_DATA_MULTI 加载多波长实验数据
    
    fprintf('=== 加载多波长实验数据 ===\n');
    fprintf('需要加载 %d 个波长的数据文件\n', length(input_params.wavelengths));
    
    exp_data_multi = [];
    
    % 逐个选择每个波长的文件
    for i = 1:length(input_params.wavelengths)
        wavelength = input_params.wavelengths(i);
        energy = input_params.energies(i);
        
        fprintf('\n选择第 %d/%d 个文件\n', i, length(input_params.wavelengths));
        fprintf('波长: %.3f nm, 能量: %.1f eV\n', wavelength, energy);
        
        [filename, pathname] = uigetfile({'*.txt;*.dat;*.csv', '数据文件 (*.txt, *.dat, *.csv)'}, ...
                                         sprintf('选择 %.3f nm 的数据文件', wavelength));
        
        if isequal(filename, 0)
            error('未选择文件');
        end
        
        filepath = fullfile(pathname, filename);
        
        % 加载数据
        data = load(filepath);
        
        exp_data.theta = data(:, 1);
        exp_data.R = data(:, 2);
        exp_data.wavelength = wavelength;
        exp_data.energy = energy;
        exp_data.energy_index = i;
        exp_data.filename = filename;
        
        exp_data_multi = [exp_data_multi; exp_data];
        
        fprintf('加载完成: %s (%d 个数据点)\n', filename, length(exp_data.theta));
    end
    
    fprintf('\n所有数据加载完成！\n');
end