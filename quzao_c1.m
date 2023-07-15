% y=Ax,现在已知y,A,求x，y为原始信号，A为字典，x为待求稀疏矩阵
% t为稀疏度
% 字典可以根据标准拉曼光谱的峰的位置进行筛选，剔除杂峰

B2 =  csvread('F:\laman\test\lib\nai.csv');
% B2 =  csvread('F:\laman\test\gauss\DDT.csv');
B2(:,1) = [];
A2=B2';
c2 = load('F:\laman\test\test\test-nai.txt');
% c2 = load('F:\laman\test\高分子处理后数据\DDT-固体.txt');
% B2 =  csvread('F:\laman\test\gauss\p_nai.csv');
% B2(:,1) = [];
% A2=B2';
% c2 = load('F:\laman\test\萘\test-nai.txt');
% y2=c2';
num= load('F:\laman\test\test\num.txt');

[M2,N2] = size(A2); %传感矩阵A为M*N矩阵
t=N2;%稀疏系数即为峰的个数
x2=OMP_c1( c2,A2,t );
if (num<0)
    for a=1:1:N2
        if (x2(a)>1)
            x2(a)=1;
        end
        if (x2(a)<0)
            x2(a)=0;
        end
    end
end
if (num>=0)
    peakp = load('F:\laman\test\test\peakp.txt');
    for a=1:1:N2
        if (x2(a)>peakp(a))
            x2(a)=peakp(a);
        end
    end
end
q2=B2';
y=q2*x2;
tt=(65:1700); 
% tt=(200:2200); 
figure (1);
subplot(2,1,1);
plot(tt,c2);
subplot(2,1,2);
plot(tt,y);
dlmwrite('F:\laman\test\萘\萘.txt',y,'delimiter', ' ');
