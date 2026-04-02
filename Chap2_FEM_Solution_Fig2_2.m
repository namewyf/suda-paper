clear; clc; close all;

%% 1. 问题参数设置
epsilon = 2^(-5);          % 振荡系数周期参数
Omega = [0, 1];            % 求解区间
f = @(x) ones(size(x));    % 源项f≡1
A_eps = @(x) 1./(2 + cos(2*pi*x/epsilon));  % 振荡扩散系数

% 解析解（公式2.3）
u_eps = @(x) x - x.^2 + (epsilon/(2*pi)).*(0.5 - x).*sin(2*pi*x/epsilon) + ...
    (epsilon^2/(4*pi^2)).*(1 - cos(2*pi*x/epsilon));

% 网格尺寸
h=2.^(-6);
x_nodes = Omega(1):h:Omega(2);  % 网格节点
x_nodes=x_nodes';
n = length(x_nodes) - 1;         % 单元数

%绘制解析解
 x_fine = linspace(Omega(1), Omega(2), 1000);  % 细网格（用于绘制光滑曲线）
 u_exact = u_eps(x_fine);        % 解析解（细网格）

 %% 3. 有限元离散与求解
    % 3.1 刚度矩阵K组装（P1有限元，分段线性基函数）
    K = sparse(n-1, n-1);  % 稀疏矩阵存储（内部节点数：n-1）
    F = zeros(n-1, 1);     % 载荷向量
     for e = 1:n  % 遍历每个单元
        % 单元节点坐标
        x1 = x_nodes(e);
        x2 = x_nodes(e+1);
        he = x2 - x1;      % 单元长度
        
        % 单元刚度矩阵（局部2x2矩阵）
        % 基函数导数：Lambda1’=1/he, Lambda2’=-1/he
        % 积分A_eps(x)在单元上的均值（因h>ε时A_eps周期振荡，用中点积分近似）
        x_mid = (x1 + x2)/2;
        A_avg = A_eps(x_mid);
        K_e = (A_avg / he) * [1, -1; -1, 1];
        
        % 单元载荷向量（积分f(x)*Lambda_i(x)，f≡1）
        F_e = (he/2) * [1; 1];  % 分段线性基函数积分结果
        
        % 组装到全局矩阵（内部节点索引：e和e+1，排除边界节点1和n+1）
        if e > 1  % 左节点为内部节点
            K(e-1, e-1) = K(e-1, e-1) + K_e(1, 1);
            F(e-1) = F(e-1) + F_e(1);
            if e < n  % 右节点为内部节点（非最后一个单元）
                K(e-1, e) = K(e-1, e) + K_e(1, 2);
            end
        end
        if e < n  % 右节点为内部节点
            K(e, e) = K(e, e) + K_e(2, 2);
            F(e) = F(e) + F_e(2);
            if e > 1  % 左节点为内部节点（非第一个单元）
                K(e, e-1) = K(e, e-1) + K_e(2, 1);
            end
        end
    end
    
    % 3.2 求解线性方程组Ku=F
    u_nodal = K \ F;  % 内部节点的数值解
    u_FEM=[0;u_nodal;0];
    plot(x_nodes,u_FEM,'b-', 'LineWidth', 2, 'DisplayName', ['数值解 h=2^{-' num2str(log2(1/h)) '}'])
hold on
     plot(x_fine, u_exact, 'r-', 'LineWidth', 2, 'DisplayName', '解析解');