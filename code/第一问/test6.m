function test6
%第一问结果：炉温曲线与数据计算
clc,clear,close all
global rol T h ms lemda
ke = 144;
ks = 10.5;
[p,p8,p9] = get_kall(ke,ks,173,198,230,257);
% [p] = get_keks(ke,ks);
% p6 = get_p6(ke,ks);
% p7 = get_p7(ke,ks);
% p8 = get_p8(ke,ks);
% p9 = get_p9(ke,ks);
f1 = @(x)p(1,1)*x^4+p(1,2)*x^3+p(1,3)*x^2+p(1,4)*x+p(1,5);%175-175
f2 = @(x)p(2,1)*x^4+p(2,2)*x^3+p(2,3)*x^2+p(2,4)*x+p(2,5);%175-195
f3 = @(x)p(3,1)*x^4+p(3,2)*x^3+p(3,3)*x^2+p(3,4)*x+p(3,5);%195-235
f4 = @(x)p(4,1)*x^4+p(4,2)*x^3+p(4,3)*x^2+p(4,4)*x+p(4,5);%235-255
f5 = @(x)p(5,1)*x^4+p(5,2)*x^3+p(5,3)*x^2+p(5,4)*x+p(5,5);%255-255
% f6 = @(x)p(6,1)*x^4+p(6,2)*x^3+p(6,3)*x^2+p(6,4)*x+p(6,5);%255-25
f6 = @(x)p6(1)*x^9+p6(2)*x^8+p6(3)*x^7+p6(4)*x^6+p6(5)*x^5+p6(6)*x^4+p6(7)*x^3+p6(8)*x^2+p6(9)*x+p6(10);
% 炉前温度
f8 = @(x)p8(1)*x^9+p8(2)*x^8+p8(3)*x^7+p8(4)*x^6+p8(5)*x^5+p8(6)*x^4+p8(7)*x^3+p8(8)*x^2+p8(9)*x+p8(10);
f9 = @(x)p9(1)*x^9+p9(2)*x^8+p9(3)*x^7+p9(4)*x^6+p9(5)*x^5+p9(6)*x^4+p9(7)*x^3+p9(8)*x^2+p9(9)*x+p9(10);

L = 50+30.5*11+5*10;
v = 78/60;%cm/s
t_real = L/v;
h = 1e-4;
ke = 32.68;
T = 0.01;
Tmax = fix(t_real/T);
lemda = 98.9;
rol = 2330;
c = get_c(25+273);%初始比热容
k = lemda/(rol*c);
m=0.15;%实际长度
ms = m*100;
mm(1) = 25;
for i = 1 : 10
    mm = [mm 30.5 5];
end
mm = [mm 30.5 25];%总距离
mm_fix = mm/v*100;
mm_real(1)=mm_fix(1);
for i = 2 : length(mm_fix)
    mm_real(i) = mm_fix(i)+mm_real(i-1);
end
mm_real = fix(mm_real);%到每层时间
mm_real(:) = mm_real(:) + 1;
r = T * k /  h^2;
A = zeros(ms+1,ms+1);
A(1,1:2) = [h*ks/lemda+1 -1];
f = zeros(1,ms+1);
for i = 2 : ms
   A(i,i-1) = -r;
   A(i,i) = 1+2*r;
   A(i,i+1) = -r;
   f(i) = 25;
end
A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
f(1) = h*ks*25/lemda;
f(ms+1) = h*ks*25/lemda;


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
    if (i>=2 && i<=mm_real(1))
        %炉前
        l_now = i/100*v;
        u0 = f8(l_now);
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ks*u0/lemda;
    elseif (i>mm_real(1) && i<=mm_real(2))
        %温区1
        u0 = 173;
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(2) && i<=mm_real(3))
        %空隙1
        l_now = (i-mm_real(2))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(3) && i<=mm_real(4))
        %温区2
        u0 = 173;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(4) && i<=mm_real(5))
        %空隙2
        l_now = (i-mm_real(4))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(5) && i<=mm_real(6))
        %温区3
        u0 = 173;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(6) && i<=mm_real(7))
        %空隙3
        l_now = (i-mm_real(6))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(7) && i<=mm_real(8))
        %温区4
        u0 = 173;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(8) && i<=mm_real(9))
        %空隙4
        l_now = (i-mm_real(8))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(9) && i<=mm_real(10))
        %温区5
        u0 = 173;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(10) && i<=mm_real(11))
        %空隙5
        l_now = (i-mm_real(10))/100*v;
        u0 = f2(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(11) && i<=mm_real(12))
        %温区6
        u0 = 198;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(12) && i<=mm_real(13))
        %空隙6
        l_now = (i-mm_real(12))/100*v;
        u0 = f3(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(13) && i<=mm_real(14))
        %温区7
        u0 = 230;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(14) && i<=mm_real(15))
        %空隙7
        l_now = (i-mm_real(14))/100*v;
        u0 = f4(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(15) && i<=mm_real(16))
        %温区8
        u0 = 257;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(16) && i<=mm_real(17))
        %空隙8
        l_now = (i-mm_real(16))/100*v;
        u0 = f5(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(17) && i<=mm_real(18))
        %温区9
        u0 = 257;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
%     elseif (i>mm_real(18) && i<=mm_real(19))
%         %空隙9
%         l_now = (i-mm_real(18))/100*v;
%         u0 = f6(l_now);
%         A(1,1:2) = [h*ke/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ke*257/lemda;
%         f(ms+1) = h*ks*u0/lemda;
%     elseif (i>mm_real(19) && i<=mm_real(20))
%         %温区10
%         u1 = 25;
% %         l_now = (i-mm_real(19))/100*70/60;
% %         u0 = f7(l_now);
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u0/lemda;
%         f(ms+1) = h*ks*u1/lemda;
%     elseif (i>mm_real(20) && i<=mm_real(21))
%         %空隙10
%         u0 = 25;
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u0/lemda;
%         f(ms+1) = h*ks*u0/lemda;
%     elseif (i>mm_real(21) && i<=mm_real(22))
%         %温区11
%         u0 = 25;
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u0/lemda;
%         f(ms+1) = h*ks*u0/lemda;
%     elseif (i>mm_real(22))
%         %炉后
%         u0 = 25;
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u0/lemda;
%         f(ms+1) = h*ks*u0/lemda;
    elseif i >mm_real(18) && i <=mm_real(21)
        l_now = (i-mm_real(18))/100*70/60;
        u0 = f9(l_now);
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ks*u0/lemda;
    elseif i >mm_real(21)
        u0 = 25;
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ks*u0/lemda;
    end
    c_now = get_c(u0+273);%当前时刻比热容
    [a,b,c] = get_abc(c_now,A);
    ytemp = chase(a,b,c,f);
    u(i,:) = ytemp;
end
q = u(:,fix((1+end)/2));
idx = find(q(:)<30);
q(idx) = [];
x = v/100:v/100:435.5+v/100;
x(idx) = [];
x = x/v;
plot(x,q,'linewidth',1)
title('炉温曲线')
xlabel('时间/s')
ylabel('温度/摄氏度')
%% 四个位置处温度
% x = x*v;
% x_3 = 25+30.5*2+5*2+30.5*0.5;%第3小温区中点
% x_6 = 25+30.5*5+5*5+30.5*0.5;%第6小温区中点
% x_7 = 25+30.5*6+5*6+30.5*0.5;%第7小温区中点
% x_8 = 25+35.5*7+30.5;%第8小温区结束处
% xx = [x_3,x_6,x_7,x_8];
% xx = xx/v;
% idx = find(abs(x-x_3)<1e-2);
% hold on
% plot(x(idx(1)),q(idx(1)),'*')
% r = q(idx)
% idx = find(abs(x-x_6)<1e-2);
% hold on
% plot(x(idx(1)),q(idx(1)),'*')
% r = q(idx)
% idx = find(abs(x-x_7)<1e-2);
% hold on
% plot(x(idx(1)),q(idx(1)),'*')
% r = q(idx)
% idx = find(abs(x-x_8)<1e-2);
% hold on
% plot(x(idx(1)),q(idx(1)),'*')
% r = q(idx)
%% 存数据
% xxx = [];
% yyy = [];
% for i = 1 : 670
%     xxx = [xxx 0.5*i];
%     idx = find(abs(x-0.5*i)<1e-2);
%     if sum(idx)>1
%        idx = idx(1) ;
%     end
%     yyy = [yyy q(idx)];
% end
% length(xxx)
% length(yyy)
% xxx = [0 xxx];
% yyy = [25 yyy];
% filename = '../result.csv';
% x = [xxx' yyy']
% % % a = csvread(filename);
% csvwrite(filename,x,1,0)
%% 三维;
% figure
% t = 0 : 1e-2 : Tmax*1e-2;
% xxx = 0 : 0.01 : ms/100;
% [xx,tt] = meshgrid(xxx,t);
% mesh(xx,tt,u)
% title('电子元件三维图像');
% xlabel('空间/ms')
% ylabel('时间/s')
% zlabel('温度/摄氏度')



function y = get_c(T)
%传入当前温度，获得当前芯片比热容
%比热容-温度多项式系数
% x = [100 200 300 400 600 800 1000 1200 1500];
% y = [259 556 712 790 867 913 946 967 992];
% p = polyfit(x,y,6);
% y = polyval(p,T);
p_t = 1e2*[-0.000000000000000   0.000000000000111  -0.000000000252512   0.000000299838736  -0.000197353916376   0.071405959681669  -2.853032199050470];
y = polyval(p_t,T);

function [a,b,c] = get_abc(cc,A)
global rol T h ms lemda
%输入比热容，获得abc
k = lemda/(rol*cc);
r = T * k /  h^2;
for i = 2 : ms
   A(i,i-1) = -r;
   A(i,i) = 1+2*r;
   A(i,i+1) = -r;
end
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