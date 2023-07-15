% y=Ax,现在已知y,A,求x，y为原始信号，A为字典，x为待求稀疏矩阵
% t为稀疏度
% 字典可以根据标准拉曼光谱的峰的位置进行筛选，剔除杂峰

B1 =  csvread('F:\laman\database_making\lib\硫.csv');
B1(:,1) = [];
A1=B1';
B2 =  csvread('F:\laman\test\test\p.csv');
B2(:,1) = [];
A2=B2';
c2 = load('F:\laman\test\test\signal.txt');
y2=c2';

[M1,N1] = size(A1); %传感矩阵A为M*N矩阵
t1=N1;%稀疏系数即为峰的个数
x1=OMP_c1( c2,A1,t1 );

for a=1:1:N1
    if (x1(a)>0.9530445840462758)
        x1(a)=0.9530445840462758;
    end
end

q1=B1';
y1=q1*x1;
tt=(70:1100); 
figure (1);
subplot(2,1,1);
plot(tt,c2);
subplot(2,1,2);
plot(tt,y1);
% dlmwrite('F:\laman\test\test\omp_signal1.txt',y1,'delimiter', ' ');

[M2,N2] = size(A2); %传感矩阵A为M*N矩阵
t2=N2;%稀疏系数即为峰的个数
x2=OMP_c1( c2,A2,t2 );

% for a=1:1:N2
%     if (x2(a)>0.9530445840462758)
%         x2(a)=0.9530445840462758;
%     end
% end

q2=B2';
y2=q2*x2;
figure (2);
subplot(2,1,1);
plot(tt,c2);
subplot(2,1,2);
plot(tt,y2);
dlmwrite('F:\laman\test\test\omp_signal3.txt',y2,'delimiter', ' ');
