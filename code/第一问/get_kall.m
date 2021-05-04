function [p,p8,p9] = get_kall(ke,ks,u_1,u_2,u_3,u_4)
%功能：传入ke、ks，输出各间层温度随位置变化拟合系数
%p——1到6情况的系数——4次拟合
%p8——炉前的系数——9次拟合
u00 = [u_1 u_1 u_2 u_3 u_4 u_4 25];
uee = [u_1 u_2 u_3 u_4 u_4 25 u_1];
% u00 = [175 175 195 235 255 255 25];
% uee = [175 195 235 255 255 25 175];
p = [];%多项式系数
val = [];%误差大小
t = [];%截断时间
p8 = [];
%% 循环遍历
for ii = 1 : 5
u0 = u00(ii);
ue = uee(ii);
h = 1e-4;
T = 0.01;
Tmax = 5000;
lemda = 0.037;
rol = 0.7833;
c = 1.02*1e3;
k = lemda/(rol*c);
r = k*T/(h^2);
ms = 500;
A = zeros(ms+1,ms+1);
A(1,1:2) = [h*ke/lemda+1 -1];

for i = 2 : ms
   A(i,i-1) = -r;
   A(i,i) = 1+2*r;
   A(i,i+1) = -r;
   f(i) = 25;
end
A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
f(1) = h*ke*u0/lemda;
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
    f(1) = h*ke*u0/lemda;
    f(ms+1) = h*ke*ue/lemda;
    ytemp = chase(a,b,c,f);
    u(i,:) = ytemp;
end
q = u(:,end);
y = u(end,:);
x = 0.01:0.01:5.01;
pp = polyfit(x,y,4);
yfit = polyval(pp,x);
p = [p;pp]; 
val = [val sum(abs(y-yfit).^2)];
end
%% p8
h = 1e-4;
T = 0.01;
Tmax = 5000;
lemda = 0.037;
rol = 0.7833;
c = 1.02*1e3;
k = lemda/(rol*c);
r = k*T/(h^2);
u0 = u00(7);
ue = uee(7);
ms = 2500;
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
y = u(end,:);
x = 0.01:0.01:25.01;
pp = polyfit(x,y,9);
yfit = polyval(pp,x);

p8 = pp;
val = sum(abs(y-yfit).^2);
%% p9
%整体 空隙9~最后
% l = 5+30.5+5+30.5+25;
l = 5+30.5+5;
h = 1e-4;
T = 0.01;
Tmax = 80000;
lemda = 0.037;
rol = 0.7833;
c = 1.02*1e3;
k = lemda/(rol*c);
r = k*T/(h^2);
u0 = u00(6);
ue = uee(6);
ms = 4050;
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
x = 0.01:0.01:40.51;
pp = polyfit(x,y,9);
yfit = polyval(pp,x);
p9 = pp;
val = sum(abs(y-yfit).^2);

% %% 前空隙
% % l = 5+30.5+5+30.5+25;
% l = 5+30.5+5;
% h = 1e-4;
% T = 0.01;
% Tmax = 80000;
% lemda = 0.037;
% rol = 0.7833;
% c = 1.02*1e3;
% k = lemda/(rol*c);
% r = k*T/(h^2);
% u0 = u_1;
% ue = u_4;
% ms = 7600;
% A = zeros(ms+1,ms+1);
% A(1,1:2) = [h*ks/lemda+1 -1];
% 
% for i = 2 : ms
%    A(i,i-1) = -r;
%    A(i,i) = 1+2*r;
%    A(i,i+1) = -r;
%    f(i) = 25;
% end
% A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
% f(1) = h*ks*u0/lemda;
% f(ms+1) = h*ks*ue/lemda;
% 
% 
% a = zeros(1,ms+1); %下对角
% b = zeros(1,ms+1); %中对角
% c = zeros(1,ms+1); %上对角
% a(1,1) = 0;
% for i = 2 : ms + 1
%     a(1,i) = A(i,i-1);
% end
% for i = 1 : ms
%     c(1,i) = A(i,i+1);
% end
% c(end) = 0;
% for i = 1 : ms + 1
%     b(1,i) = A(i,i); 
% end
% u = zeros(Tmax+1,ms+1);
% ytemp = chase(a,b,c,f);
% u(1,:) = 25;
% u(2,:) = ytemp;
% for i = 2 : Tmax + 1 
%     f = ytemp';
%     f(1) = h*ks*u0/lemda;
%     f(ms+1) = h*ks*ue/lemda;
%     ytemp = chase(a,b,c,f);
%     u(i,:) = ytemp;
% end
% q = u(:,end);
% y = u(end,:);
% x = 0.01:0.01:76.01;
% pp = polyfit(x,y,9);
% yfit = polyval(pp,x);
% p10 = pp;
% val = sum(abs(y-yfit).^2);