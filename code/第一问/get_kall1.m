function [p10] = get_kall1(ke,ks,u_1,u_2,u_3,u_4)
%功能：传入ke、ks，输出各间层温度随位置变化拟合系数
%% 前空隙
% l = 5+30.5+5+30.5+25;
l = 5+30.5+5;
h = 1e-4;
T = 0.01;
Tmax = 180000;
lemda = 0.037;
rol = 0.7833;
c = 1.02*1e3;
k = lemda/(rol*c);
r = k*T/(h^2);
u0 = u_1;
ue = u_4;
ms = 7600;
A = zeros(ms+1,ms+1);
A(1,1:2) = [h*ks/lemda+1 -1];

for i = 2 : ms
   A(i,i-1) = -r;
   A(i,i) = 1+2*r;
   A(i,i+1) = -r;
   f(i) = 25;
end
A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
f(1) = h*ks*u0/lemda;
f(ms+1) = h*ks*ue/lemda;


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
    f(ms+1) = h*ks*ue/lemda;
    ytemp = chase(a,b,c,f);
    u(i,:) = ytemp;
end
q = u(:,end);
y = u(end,:);
x = 0.01:0.01:76.01;
pp = polyfit(x,y,9);
figure
plot(y)
yfit = polyval(pp,x);
p10 = pp;