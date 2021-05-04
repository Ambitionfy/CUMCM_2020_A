function test_PSO4
clear,clc,close all
%粒子群算法
s1 = 1.49445;
s2 = 1.49445;%学习因子
Maxg = 60;%进化次数
SizePop = 100;%种群规模
x1max = 185;
x1min = 165;
x2max = 205;
x2min = 185;
x3max = 245;
x3min = 225;
x4max = 265;
x4min = 245;
x5max = 100/60;
x5min = 65/60;
Vmax = 1;
Vmin = -1;
% Xmax = 1200;
% Xmin = 800;
% Ymax = 1000;
% Ymin = 700;
par_num = 5;
%产生初始粒子和速度
for i = 1 : SizePop
   Pop(i,1) = 175+10*rands(1); 
   Pop(i,2) = 195+10*rands(1);
   Pop(i,3) = 235+10*rands(1); 
   Pop(i,4) = 255+10*rands(1);
   Pop(i,5) = 82.5/60+17.5/60*rands(1); 

   V(i,:) = rands(1,par_num);
   %计算适应度
   
   Fitness(i) = fun1(Pop(i,1),Pop(i,2),Pop(i,3),Pop(i,4),Pop(i,5));
end

%%最优位置初始化
[BestFitness,BestIndex] = min(Fitness);
gBest = Pop(BestIndex,:);%全局最优
pBest = Pop;
FitnesspBest = Fitness;%个体最佳适应度
FitnessgBest = BestFitness;%全局最佳适应度

%%迭代寻优
for i  = 1 : Maxg
    for j = 1 : SizePop
        %种群中各粒子速度更新
        V(j,:) = V(j,:) + s1*rand*(pBest(j,:)-Pop(j,:)) + s2*rand*(gBest-Pop(j,:));
        %更新速度——超过阈值的设为阈值
        V(j,find(V(j,:) > Vmax)) = Vmax;
        V(j,find(V(j,:) < Vmin)) = Vmin;
        %种群粒子位置更新,0.5可换成0-1任意数
        Pop(j,:) = Pop(j,:) + 0.5*V(j,:);
        Pop(j,find(Pop(j,1) > x1max)) = x1max;
        Pop(j,find(Pop(j,1) < x1min)) = x1min;
        Pop(j,find(Pop(j,2) > x2max)) = x2max;
        Pop(j,find(Pop(j,2) < x2min)) = x2min;
        Pop(j,find(Pop(j,3) > x3max)) = x3max;
        Pop(j,find(Pop(j,3) < x3min)) = x3min;
        Pop(j,find(Pop(j,4) > x4max)) = x4max;
        Pop(j,find(Pop(j,4) < x4min)) = x4min;
        Pop(j,find(Pop(j,5) > x5max)) = x5max;
        Pop(j,find(Pop(j,5) < x5min)) = x5min;
        Fitness(j) = fun1(Pop(j,1),Pop(j,2),Pop(j,3),Pop(j,4),Pop(j,5));%计算适应度
        %个体最优更新
        if Fitness(j) < FitnesspBest(j)
           pBest(j,:) =  Pop(j,:);
           FitnesspBest(j) = Fitness(j);
        end
        %群体最优更新
        if Fitness(j) < FitnessgBest
           gBest = Pop(j,:);
           FitnessgBest = Fitness(j);
        end
    end
    yy(i) = FitnessgBest;
end
%%输出结果
[FitnessgBest, gBest]
% [x,y] = meshgrid(800:1200,700:1000);
% z = fun1(x,y);
% plot3(x,y,z);
% hold on 
% plot3(gBest(1), gBest(2), FitnessgBest, 'bo','linewidth', 1.5)
%  
figure
plot(yy)
title('最优个体适应度', 'fontsize', 12);
xlabel('进化代数', 'fontsize', 12);
ylabel('适应度', 'fontsize', 12);

