function [ x ] = OMP_c1( y,A,t )
    % y=Ax,现在一直y,A,求x
    % t为稀疏度
    [M,N] = size(A); %传感矩阵A为M*N矩阵
    x = zeros(N,1); %用来存储恢复的x(列向量)
    At = zeros(M,t); %用来存储迭代过程中A被选择的列
    Pos_x = zeros(1,t); %用来存储迭代过程中A被选择的列序号
    r_n = y; %初始化残差为y
    for ii=1:t %迭代t次
        product = A'*r_n;% A各列与残差的内积
        [val,pos] = max(abs(product));%找到最大内积绝对值
        At(:,ii) = A(:,pos);%存储这一列
        Pos_x(ii) = pos;%存储这一列的序号
        A(:,pos) = zeros(M,1);%清零A的这一列，因为它与残差正交
        %y=At(:,1:ii)*x，以下求x的最小二乘解
        x_ls = (At(:,1:ii)'*At(:,1:ii))^(-1)*At(:,1:ii)'*y;%最小二乘解
        %At(:,1:ii)*x_ls是y在At(:,1:ii)列空间上的正交投影
        r_n = y - At(:,1:ii)*x_ls;%更新残差        
    end
    x(Pos_x)=x_ls;%恢复出的x
end
