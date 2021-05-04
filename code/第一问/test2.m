function test2
clear,clc,close all
%改变u0，ue，获得间隙值
h = 1e-4;
ke=144;   
ks = 10.5;
T = 0.01;
Tmax = 2000;
lemda = 0.037;
rol = 0.7833;
c = 1.02*1e3;
k = lemda/(rol*c);
r = k*T/(h^2);
u0 = 255;
ue = 255;
c1 = [-lemda lemda+lemda -lemda]/h; %1间隙边界条件
ms = 500;
A = zeros(ms+1,ms+1);
A(1,1:2) = [h*ks/lemda+1 -1];

for i = 2 : ms
   A(i,i-1) = -r;
   A(i,i) = 1+2*r;
   A(i,i+1) = -r;
   f(i) = 25;
end
A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
f(1) = h*ks*u0/lemda;
f(ms+1) = h*ke*ue/lemda;


a = zeros(1,ms+1); %下对角
b = zeros(1,ms+1); %中对角
c = zeros(1,ms+1); %上对角
a(1,1) = 0;
for i = 2 : ms + 1
    a(1,i) = A(i,i-1);
end
for i = 1 : ms
    c(1,i) = A(i,i+1);
end
c(end) = 0;
for i = 1 : ms + 1
    b(1,i) = A(i,i); 
end
u = zeros(Tmax+1,ms+1);
ytemp = chase(a,b,c,f);
u(1,:) = 25;
u(2,:) = ytemp;
for i = 2 : Tmax + 1 
    f = ytemp';
    f(1) = h*ks*u0/lemda;
    f(ms+1) = h*ke*ue/lemda;
    ytemp = chase(a,b,c,f);
    u(i,:) = ytemp;
end
q = u(:,end);
t = 0 : Tmax;
xxx = 0 : ms;
[xx,tt] = meshgrid(xxx,t);
figure(1)
mesh(xx,tt,u)
title('温度随时间、空间分布图')
xlabel('空间/0.1mm')
ylabel('时间/s')
zlabel('温度/摄氏度')
t = 0;
for i = 2 : length(q) 
   if abs(q(i)-q(i-1))<1e-3
       t = i;
       break
   end
end
t
figure(2)
plot(0.01:0.01:length(q)/100,q,'linewidth',1)
title('边界温度随时间变化图')
xlabel('时间/s')
ylabel('温度/摄氏度')
figure(3)
plot(0.01:0.01:5.01,u(t,:),'linewidth',1)
title('短时间后空间随温度变化曲线')
xlabel('空间/cm')
ylabel('温度/摄氏度')
y = u(t,:);
x = 0.01:0.01:5.01;
p = polyfit(x,y,9);
yfit = polyval(p,x);
figure(4)
plot(x,y,'b-.',x,yfit,'r')
title('拟合曲线')
xlabel('位置/cm')
ylabel('摄氏度')
legend('原始数据','拟合数据')
p
val = sum(abs(y-yfit).^2)
figure
plot(u(:,fix((1+end)/2)))
