clear; clc; close all;

%% 问题参数设置
%HINT 还需要在单刚矩阵部分修改A_eps_inv_t
epsilon = 2^(-5);          % 振荡系数周期参数
Omega = [0, 1];            % 求解区间
f = @(x) ones(size(x));    % 源项f≡1
% A_eps = @(x) 1./(2 + cos(2*pi*x/epsilon));  % 振荡扩散系数
% A_eps_inv = @(x) 1 ./ A_eps(x); % A_eps的逆
A_eps = @(x) 1./(2 -sin(2*pi*tan(15*pi*x/32)));  % 扩散系数(2.17)，
A_eps_inv = @(x) 1 ./ A_eps(x); % A_eps的逆



% 网格尺寸
H=2.^(-2);
h=2.^(-6);
x_nodes = Omega(1):H:Omega(2);  % 网格节点
x_nodes=x_nodes';
n = length(x_nodes) - 1;         % 单元数
ratio=H/h;
LeftBasisValue=zeros((ratio-1)*n,1);
RightBasisValue=zeros((ratio-1)*n,1);
InteriorNodes=zeros((ratio-1)*n,1);



 %% LOD有限元离散与求解
    % 3.1 刚度矩阵K组装（Λ_j^{ms}）
    K = sparse(n-1, n-1);  % 稀疏矩阵存储（内部节点数：n-1）
    F = zeros(n-1, 1);     % 载荷向量
     for e = 1:n  % 遍历每个单元
        % 单元节点坐标
        x1 = x_nodes(e);
        x2 = x_nodes(e+1);
        he = x2 - x1;      % 单元长度
        
        % 单元刚度矩阵（局部2x2矩阵）
    
        c=integral(A_eps_inv, x1, x2);
        K_e = (1 / c) * [1, -1; -1, 1];
        
        % 单元载荷向量（积分f(x)*Λ_j^{ms}(x)，f≡1）
      % 用符号运算，后面可替换成数值计算
       syms x t
       A_eps_inv_t=2 -sin(2*pi*tan(15*pi*t/32));% 扩散系数(2.17)的逆
      % A_eps_inv_t=2 + cos(2*pi*t/epsilon); % 振荡扩散系数(2.2)的逆
       f1=int(A_eps_inv_t,t,x1,x);
       F_e =[int(1-f1/c,x, x1, x2);int(f1/c, x,x1, x2)];
       

       %细网格上离散点基函数的值

       LeftBasisValue((ratio-1)*(e-1)+1:(ratio-1)*e)=1-double(subs(f1,x,x1+h:h:x2-h)/c);
       InteriorNodes((ratio-1)*(e-1)+1:(ratio-1)*e)=x1+h:h:x2-h;
       RightBasisValue((ratio-1)*(e-1)+1:(ratio-1)*e)=double(subs(f1,x,x1+h:h:x2-h)/c);


        
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

     %  plot 基函数Λ_j^{ms}(x)
     figure(1)
pointidx=3;
idx=(pointidx-2)*(ratio-1)+1:(pointidx-1)*(ratio-1);
RightPart= RightBasisValue(idx);
idx=(pointidx-1)*(ratio-1)+1:(pointidx)*(ratio-1);
LeftPart= LeftBasisValue(idx);
xelem=(pointidx-2)*H:h:(pointidx)*H;
basiselem=[0;RightPart;1;LeftPart;0];
plot(xelem,basiselem,'b-')
axis([0 1 0 1])
    
    % 3.2 求解线性方程组Ku=F
    u_nodal = K \ F;  % TH内部节点的数值解
    u_FEM=[0;u_nodal;0];
    Left=[0;u_nodal];
    Left=repmat(Left,1,ratio-1);
    Left=Left';
    Left=Left(:);
    Right=[u_nodal;0];
    Right=repmat(Right,1,ratio-1);
    Right=Right';
    Right=Right(:);
    InteriorNodesValue=LeftBasisValue.*Left+RightBasisValue.*Right;
    x_nodes_h=[x_nodes;InteriorNodes];
    u_FEM_ms=[u_FEM;InteriorNodesValue];%TH每个单元内部Th离散点值

    [x_nodes_h,idx]=sort(x_nodes_h);% Th节点按顺序排列
    u_FEM_ms=u_FEM_ms(idx);% Th节点按顺序排列
figure(2)
    plot(x_nodes_h,u_FEM_ms,'b-', 'LineWidth', 2, 'DisplayName', ['数值解 H=2^{-' num2str(log2(1/H)) '}'])
hold on 
%% 细网格上有限元离散与求解
x_nodes = Omega(1):h:Omega(2);  % 网格节点
x_nodes=x_nodes';
n = length(x_nodes) - 1;         % 单元数
    % 刚度矩阵K组装（P1有限元，分段线性基函数）
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
    plot(x_nodes,u_FEM,'r-*', 'LineWidth', 2, 'DisplayName', ['数值解 h=2^{-' num2str(log2(1/h)) '}'])
  
