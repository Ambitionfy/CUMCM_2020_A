function test_PSO4
clear,clc,close all
%����Ⱥ�㷨
s1 = 1.49445;
s2 = 1.49445;%ѧϰ����
Maxg = 60;%��������
SizePop = 100;%��Ⱥ��ģ
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
%������ʼ���Ӻ��ٶ�
for i = 1 : SizePop
   Pop(i,1) = 175+10*rands(1); 
   Pop(i,2) = 195+10*rands(1);
   Pop(i,3) = 235+10*rands(1); 
   Pop(i,4) = 255+10*rands(1);
   Pop(i,5) = 82.5/60+17.5/60*rands(1); 

   V(i,:) = rands(1,par_num);
   %������Ӧ��
   
   Fitness(i) = fun1(Pop(i,1),Pop(i,2),Pop(i,3),Pop(i,4),Pop(i,5));
end

%%����λ�ó�ʼ��
[BestFitness,BestIndex] = min(Fitness);
gBest = Pop(BestIndex,:);%ȫ������
pBest = Pop;
FitnesspBest = Fitness;%���������Ӧ��
FitnessgBest = BestFitness;%ȫ�������Ӧ��

%%����Ѱ��
for i  = 1 : Maxg
    for j = 1 : SizePop
        %��Ⱥ�и������ٶȸ���
        V(j,:) = V(j,:) + s1*rand*(pBest(j,:)-Pop(j,:)) + s2*rand*(gBest-Pop(j,:));
        %�����ٶȡ���������ֵ����Ϊ��ֵ
        V(j,find(V(j,:) > Vmax)) = Vmax;
        V(j,find(V(j,:) < Vmin)) = Vmin;
        %��Ⱥ����λ�ø���,0.5�ɻ���0-1������
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
        Fitness(j) = fun1(Pop(j,1),Pop(j,2),Pop(j,3),Pop(j,4),Pop(j,5));%������Ӧ��
        %�������Ÿ���
        if Fitness(j) < FitnesspBest(j)
           pBest(j,:) =  Pop(j,:);
           FitnesspBest(j) = Fitness(j);
        end
        %Ⱥ�����Ÿ���
        if Fitness(j) < FitnessgBest
           gBest = Pop(j,:);
           FitnessgBest = Fitness(j);
        end
    end
    yy(i) = FitnessgBest;
end
%%������
[FitnessgBest, gBest]
% [x,y] = meshgrid(800:1200,700:1000);
% z = fun1(x,y);
% plot3(x,y,z);
% hold on 
% plot3(gBest(1), gBest(2), FitnessgBest, 'bo','linewidth', 1.5)
%  
figure
plot(yy)
title('���Ÿ�����Ӧ��', 'fontsize', 12);
xlabel('��������', 'fontsize', 12);
ylabel('��Ӧ��', 'fontsize', 12);

