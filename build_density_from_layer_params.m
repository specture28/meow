function grid_model = build_density_from_layer_params(layer_params, input_params)
% BUILD_DENSITY_FROM_LAYER_PARAMS 全局连续性方法构建密度分布

    % 固定网格设置
    fixed_depth = 15.5;
    z_grid = 0:input_params.dz:fixed_depth;
    num_points = length(z_grid);
    
    % 初始化所有元素密度
    all_elements = {'Si', 'Co', 'Mg', 'O', 'C'};
    density_grid = struct();
    for i = 1:length(all_elements)
        density_grid.(all_elements{i}) = zeros(size(z_grid));
    end
    
    % 解析层参数
    param_idx = 1;
    layer_info = [];  % 存储所有层的信息
    
    % 第一步：解析所有层参数，构建层结构
    for layer = 1:input_params.num_layers
        layer_elements = input_params.layer_info{layer}.elements;
        layer_thickness = input_params.layer_info{layer}.thickness;
        
        if ismember('SiO2', layer_elements)
            % === SiO2层：分解为Si+O ===
            
            % 🔧 修正：按照construct_initial_layer_params.m的参数顺序解析
            rho_Si = layer_params(param_idx);           % 1. Si密度
            sigma_L_shared = layer_params(param_idx+1); % 2. 共享σL  
            sigma_R_Si = layer_params(param_idx+2);     % 3. Si独立σR
            z0_shared = layer_params(param_idx+3);      % 4. 共享z0
            thickness_shared = layer_params(param_idx+4); % 5. 共享thickness
            rho_O = layer_params(param_idx+5);          % 6. O密度
            sigma_R_O = layer_params(param_idx+6);      % 7. O独立σR
            param_idx = param_idx + 7;
            
            % 计算层边界
            layer_start = z0_shared - thickness_shared/2;
            layer_end = z0_shared + thickness_shared/2;
            
            % Si组分
            si_layer = struct();
            si_layer.element = 'Si';
            si_layer.rho = rho_Si;
            si_layer.sigma_L = sigma_L_shared;
            si_layer.sigma_R = sigma_R_Si;  % Si独立右界面
            si_layer.z0 = z0_shared;
            si_layer.thickness = thickness_shared;
            si_layer.layer_start = layer_start;
            si_layer.layer_end = layer_end;
            layer_info = [layer_info, si_layer];
            
            % O组分  
            o_layer = struct();
            o_layer.element = 'O';
            o_layer.rho = rho_O;
            o_layer.sigma_L = sigma_L_shared;  % 共享左界面
            o_layer.sigma_R = sigma_R_O;       % O独立右界面
            o_layer.z0 = z0_shared;
            o_layer.thickness = thickness_shared;
            o_layer.layer_start = layer_start;
            o_layer.layer_end = layer_end;
            layer_info = [layer_info, o_layer];
            
            % 处理SiO2层中的其他元素
            other_elements = setdiff(layer_elements, {'SiO2'});
            for elem_idx = 1:length(other_elements)
                elem = other_elements{elem_idx};
                
                rho = layer_params(param_idx);
                sigma_L = layer_params(param_idx+1);
                sigma_R = layer_params(param_idx+2);
                z0 = layer_params(param_idx+3);
                param_idx = param_idx + 4;
                
                elem_layer = struct();
                elem_layer.element = elem;
                elem_layer.rho = rho;
                elem_layer.sigma_L = sigma_L;
                elem_layer.sigma_R = sigma_R;
                elem_layer.z0 = z0;
                elem_layer.thickness = layer_thickness;  % 使用原始层厚度
                elem_layer.layer_start = z0 - layer_thickness/2;
                elem_layer.layer_end = z0 + layer_thickness/2;
                layer_info = [layer_info, elem_layer];
            end
            
        else
            % === 普通层 ===
            for elem_idx = 1:length(layer_elements)
                elem = layer_elements{elem_idx};
                
                rho = layer_params(param_idx);
                sigma_L = layer_params(param_idx+1);
                sigma_R = layer_params(param_idx+2);
                z0 = layer_params(param_idx+3);
                param_idx = param_idx + 4;
                
                elem_layer = struct();
                elem_layer.element = elem;
                elem_layer.rho = rho;
                elem_layer.sigma_L = sigma_L;
                elem_layer.sigma_R = sigma_R;
                elem_layer.z0 = z0;
                elem_layer.thickness = layer_thickness;
                elem_layer.layer_start = z0 - layer_thickness/2;
                elem_layer.layer_end = z0 + layer_thickness/2;
                layer_info = [layer_info, elem_layer];
            end
        end
    end
    
    % 提取基底参数
    total_film_params = param_idx - 1;
    if length(layer_params) > total_film_params
        substrate_sigma_L = layer_params(total_film_params + 1);
    else
        substrate_sigma_L = input_params.substrate_sigma_L;
    end
    
    % 找到最深的层边界作为基底界面
    substrate_interface_pos = 0;
    for i = 1:length(layer_info)
        substrate_interface_pos = max(substrate_interface_pos, layer_info(i).layer_end);
    end
    
    % 第二步：对每个深度点，计算所有元素的密度
    for z_idx = 1:num_points
        z = z_grid(z_idx);
        
        % 对每个元素，累加所有相关层的贡献
        for elem_name = all_elements
            elem = elem_name{1};
            total_density = 0;
            
            % 累加所有包含该元素的层贡献
            for layer_idx = 1:length(layer_info)
                layer_data = layer_info(layer_idx);
                if strcmp(layer_data.element, elem)
                    
                    % 计算该层对当前深度的贡献
                    erf_start = erf((z - layer_data.layer_start) / (layer_data.sigma_L * sqrt(2)));
                    erf_end = erf((z - layer_data.layer_end) / (layer_data.sigma_R * sqrt(2)));
                    
                    contribution = layer_data.rho/2 * (erf_start - erf_end);
                    total_density = total_density + max(0, contribution);
                end
            end
            
            % 添加基底贡献
            if strcmp(elem, 'Si')
                % Si基底：深层密度趋向substrate_density
                erf_substrate = erf((z - substrate_interface_pos) / (substrate_sigma_L * sqrt(2)));
                substrate_contribution = input_params.substrate_density * (1 + erf_substrate) / 2;
                total_density = total_density + substrate_contribution;
                
            elseif strcmp(elem, 'O')
                % O基底：深层密度趋向0（自然处理，无需额外代码）
                % 因为基底不含O，所以不添加任何贡献
                
            else
                % 其他元素：深层密度趋向0（自然处理）
                % 因为基底不含这些元素，所以不添加任何贡献
            end
            
            density_grid.(elem)(z_idx) = total_density;
        end
    end
    
    % 构建输出结构
    grid_model = struct();
    grid_model.z_grid = z_grid;
    grid_model.density_grid = density_grid;
    grid_model.num_points = num_points;
    grid_model.dz = input_params.dz;
    grid_model.all_elements = all_elements;
    grid_model.calculation_depth = fixed_depth;
    
    grid_model.substrate = struct();
    grid_model.substrate.material = 'Si';
    grid_model.substrate.density = input_params.substrate_density;
    grid_model.substrate.interface_position = substrate_interface_pos;
    grid_model.substrate.sigma_L = substrate_sigma_L;
    
    grid_model.layer_info_decomposed = layer_info;  % 保存分解后的层信息
    grid_model.sio2_decomposed = true;
    
    fprintf('✅ 密度构建完成：%d个分解层，基底界面位置=%.2f nm\n', ...
            length(layer_info), substrate_interface_pos);
    
    % 验证边界条件
    surface_si = density_grid.Si(1);
    surface_o = density_grid.O(1);
    deep_si = density_grid.Si(end);
    deep_o = density_grid.O(end);
    
    fprintf('📊 边界检查：表面Si=%.3f O=%.3f, 深层Si=%.3f O=%.3f\n', ...
            surface_si, surface_o, deep_si, deep_o);
end