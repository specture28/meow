function R = parratt(theta, lambda, layers)
% PARRATT 使用Parratt递归法计算X射线反射率
%   R = PARRATT(theta, lambda, layers) 计算给定角度和波长下的反射率
%   
%   约定：
%   ñ = 1 - δ - iβ  (X射线光学标准约定)
%   n = sqrt(1 - χ)  # 对应此约定的χ定义
%
%   输入参数:
%   theta - 入射角度数组(度)
%   lambda - X射线波长(nm)
%   layers - 包含层信息的结构体数组，每层包含:
%            .thickness - 厚度(nm)
%            .chi - 极化率(复数)
%
%   输出:
%   R - 反射率数组

    % 检查输入参数
    if nargin < 3
        error('需要提供角度、波长和层结构');
    end
    
    % 转换角度为弧度
    theta_rad = theta * pi / 180;
    
    % 波矢
    k0 = 2*pi/lambda;
    
    % 层数
    N = length(layers);
    
    % 结果数组
    R = zeros(size(theta_rad));
    
    % 处理每个角度
    for i = 1:length(theta_rad)
        % 计算每层的kz
        kz = zeros(N, 1);
        for j = 1:N
            % 修正：使用标准X射线光学约定
            n = sqrt(1 - layers(j).chi);  % 改为减号
            
            % kz计算保持不变
            kz_temp = k0 * sqrt(n^2 - cos(theta_rad(i))^2);
            
            % 确保选择物理正确的分支：Im(kz) > 0
            if imag(kz_temp) < 0
                kz(j) = -kz_temp;
            else
                kz(j) = kz_temp;
            end
        end
        
        % Parratt递归
        r = 0; % 底层反射率初始化为0
        
        for j = N-1:-1:1
            % 计算菲涅尔系数
            rj = (kz(j) - kz(j+1)) / (kz(j) + kz(j+1));
            
            % 相位因子
            if j < N  % 避免对底层计算相位
                % 相位计算 - 检查是否需要符号调整
                phase_arg = kz(j+1) * layers(j+1).thickness;
                
                % 在标准X射线约定下，相位因子保持不变
                phase = exp(2i * phase_arg);
                
                % 递归公式
                r = (rj + r * phase) / (1 + rj * r * phase);
            else
                r = rj;
            end
        end
        
        % 计算反射率
        R(i) = abs(r)^2;
    end
end