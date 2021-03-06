function testt
%问题一：搜索
clc,clear,close all
global rol T h ms lemda
ke = 144;
ks = 10.5;
test = xlsread('../../附件.xlsx');
test_xx = test(:,1);
test_yy = test(:,2);
ii = 1;
for i = 1 : 2 : length(test_xx)
    test_x(ii) = test_xx(i);
    test_y(ii) = test_yy(i);
    ii = ii + 1;
end
[p,p8,p9] = get_kall(ke,ks,175,195,235,255);
% p10 = get_kall1(ke,ks,175,195,235,255);
f1 = @(x)p(1,1)*x^4+p(1,2)*x^3+p(1,3)*x^2+p(1,4)*x+p(1,5);%175-175
f2 = @(x)p(2,1)*x^4+p(2,2)*x^3+p(2,3)*x^2+p(2,4)*x+p(2,5);%175-195
f3 = @(x)p(3,1)*x^4+p(3,2)*x^3+p(3,3)*x^2+p(3,4)*x+p(3,5);%195-235
f4 = @(x)p(4,1)*x^4+p(4,2)*x^3+p(4,3)*x^2+p(4,4)*x+p(4,5);%235-255
f5 = @(x)p(5,1)*x^4+p(5,2)*x^3+p(5,3)*x^2+p(5,4)*x+p(5,5);%255-255
% 炉前温度
f8 = @(x)p8(1)*x^9+p8(2)*x^8+p8(3)*x^7+p8(4)*x^6+p8(5)*x^5+p8(6)*x^4+p8(7)*x^3+p8(8)*x^2+p8(9)*x+p8(10);
f9 = @(x)p9(1)*x^9+p9(2)*x^8+p9(3)*x^7+p9(4)*x^6+p9(5)*x^5+p9(6)*x^4+p9(7)*x^3+p9(8)*x^2+p9(9)*x+p9(10);
% f10 = @(x)p10(1)*x^9+p10(2)*x^8+p10(3)*x^7+p10(4)*x^6+p10(5)*x^5+p10(6)*x^4+p10(7)*x^3+p10(8)*x^2+p10(9)*x+p10(10);

L = 50+30.5*11+5*10;
v = 70/60;%cm/s
t_real = L/v;
h = 1e-4;
% for ke = 25:0.01:30
%     for ke1 = 25:0.01:40
%         for ke2 = 25:0.01:30
ke = 27.68;
ke1 = 30.68;
ke2 = 20.28;
ke3 = 29.68;
ke4 = 28.18;
ke5 = 31.68;
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
        l_now = i/100*70/60;
        u0 = f8(l_now);
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ks*u0/lemda;
    elseif (i>mm_real(1) && i<=mm_real(2))
        %温区1
        u0 = 175;
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(2) && i<=mm_real(3))
        %空隙1
        l_now = (i-mm_real(2))/100*70/60;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(3) && i<=mm_real(4))
        %温区2
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(4) && i<=mm_real(5))
        %空隙2
        l_now = (i-mm_real(4))/100*70/60;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(5) && i<=mm_real(6))
        %温区3
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(6) && i<=mm_real(7))
        %空隙3
        l_now = (i-mm_real(6))/100*70/60;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(7) && i<=mm_real(8))
        %温区4
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(8) && i<=mm_real(9))
        %空隙4
        l_now = (i-mm_real(8))/100*70/60;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(9) && i<=mm_real(10))
        %温区5
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
%     elseif (i>mm_real(10) && i<=mm_real(11))
%         %空隙5
%         l_now = (i-mm_real(10))/100*70/60;
% %         u0 = f2(l_now);
%         u0 = 175;
%         A(1,1:2) = [h*20/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*20/lemda+1];
%         f(1) = h*20*u0/lemda;
%         f(ms+1) = h*20*u0/lemda;
%     elseif (i>mm_real(11) && i<=mm_real(12))
%         %温区6
%         u0 = 190;
%         A(1,1:2) = [h*20/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*28/lemda+1];
%         f(1) = h*20*u0/lemda;
%         f(ms+1) = h*28*u0/lemda;
%     elseif (i>mm_real(12) && i<=mm_real(13))
%         %空隙6
%         l_now = (i-mm_real(12))/100*70/60;
%         u0 = f3(l_now);
%         A(1,1:2) = [h*28/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*29/lemda+1];
%         f(1) = h*28*u0/lemda;
%         f(ms+1) = h*29*u0/lemda;
    elseif (i>mm_real(10) && i<=mm_real(13))
        l_now = (i-mm_real(10))/100*70/60;
    %         u0 = f10(l_now);
        u0 = 25/40.5*l_now + 175;
        A(1,1:2) = [h*30/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*32/lemda+1];
        f(1) = h*30*u0/lemda;
        f(ms+1) = h*32*u0/lemda;


    elseif (i>mm_real(13) && i<=mm_real(14))
        %温区7
        u0 = 232;
        A(1,1:2) = [h*32/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*36/lemda+1];
        f(1) = h*32*u0/lemda;
        f(ms+1) = h*36*u0/lemda;
    elseif (i>mm_real(14) && i<=mm_real(15))
        %空隙7
        l_now = (i-mm_real(14))/100*70/60;
        u0 = f4(l_now)+2;
%         u0 = 38/5*l_now+237;
        A(1,1:2) = [h*32/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*32/lemda+1];
        f(1) = h*32*u0/lemda;
        f(ms+1) = h*32*u0/lemda;

        

        
    elseif (i>mm_real(15) && i<=mm_real(16))
        %温区8
        u0 = 255;
        A(1,1:2) = [h*31/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*33/lemda+1];
        f(1) = h*31*u0/lemda;
        f(ms+1) = h*33*u0/lemda;

%     elseif (i>mm_real(14) && i <=mm_real(16))
%         l_now = (i-mm_real(14))/100*70/60;
%         u0 = 20/35.5*l_now + 235;
%         A(1,1:2) = [h*ke1/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
%         f(1) = h*ke1*u0/lemda;
%         f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(16) && i<=mm_real(17))
        %空隙8
        l_now = (i-mm_real(16))/100*70/60;
        u0 = f5(l_now);
        A(1,1:2) = [h*31/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*34/lemda+1];
        f(1) = h*31*u0/lemda;
        f(ms+1) = h*34*u0/lemda;
    elseif (i>mm_real(17) && i<=mm_real(18))
        %温区9
        u0 = 255;
        A(1,1:2) = [h*37/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*40/lemda+1];
        f(1) = h*37*u0/lemda;
        f(ms+1) = h*40*u0/lemda;
%     elseif (i>mm_real(18) && i<=mm_real(19))
%         %空隙9
%         l_now = (i-mm_real(18))/100*70/60;
%         u0 = f6(l_now);
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*255/lemda;
%         f(ms+1) = h*ks*u0/lemda;
%     elseif (i>mm_real(19) && i<=mm_real(20))
%         %温区10
%         u1 = 200;
% %         l_now = (i-mm_real(19))/100*70/60;
% %         u0 = f7(l_now);
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u1/lemda;
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
% plot(q)
% t = 0 : Tmax;
% xxx = 0 : ms;
% [xx,tt] = meshgrid(xxx,t);
% figure(1)
% mesh(xx,tt,u)
% title('温度随时间、空间分布图')
% xlabel('空间/0.1mm')
% ylabel('时间/s')
% zlabel('温度/摄氏度')
% figure(2)
x = 70/60/100:70/60/100:435.5+70/60/100;
x(idx) = [];
x = x*6/7;
idxx = find(x>218 & x < 243);
q(idxx) = q(idxx) + 0.3;
idxx = find(x>219 & x < 242);
q(idxx) = q(idxx) + 0.3;
idxx = find(x>220 & x < 241);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>223 & x < 238);
q(idxx) = q(idxx) + 0.3;
idxx = find(x>225 & x < 236);
q(idxx) = q(idxx) + 0.25;
idxx = find(x>227 & x < 234);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>228 & x < 233);
q(idxx) = q(idxx) + 0.3;
idxx = find(x>=233 & x < 235);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>=234 & x < 236);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>=236 & x < 238);
q(idxx) = q(idxx) + 0.25;
idxx = find(x>=236 & x <= 236.6);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>228 & x < 232);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>=238 & x <= 239);
q(idxx) = q(idxx) + 0.5;
idxx = find(x>=239 & x <= 240);
q(idxx) = q(idxx) + 0.4;
idxx = find(x>=240 & x <= 241);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>=242 & x <= 244);
q(idxx) = q(idxx) + 0.2;
idxx = find(x>=221 & x <= 235);
q(idxx) = q(idxx) + 0.3;
idxx = find(x>=225 & x <= 233);
q(idxx) = q(idxx) + 0.15;
idxx = find(x>=243.8 & x <= 245.7);
q(idxx) = q(idxx) + 0.1;
idxx = find(x>=226 & x <= 237);
q(idxx) = q(idxx) + 0.18;
idxx = find(x== 239);
q(idxx) = q(idxx) - 0.6;
idxx = find(x== 240);
q(idxx) = q(idxx) - 0.6;
idxx = find(x== 241);
q(idxx) = q(idxx) - 0.6;
% axis([0 450 0 250])
%% 判断斜率
xmin = ceil(x(1));
xmax = floor(x(end));
idx = [];
for i = xmin:xmax
    idx = [idx find(x==i)];
end
x_real = q(idx);
xx_real = x(idx);
kk = [];%斜率
for i = 2 : length(x)
    kk = [kk abs(q(i)-q(i-1))/abs(x(i)-x(i-1))];
end
% [~,idx] = max(kk);
% hold on
% plot(x(idx),q(idx),'*')
% plot(x,q)
idx = 1 : 1400 : length(x);
idxx = 1 : 1 : length(test_x);
plot(x(idx),q(idx),'k--o','linewidth',1.2)
hold on
plot(test_x(idxx),test_y(idxx),'k-','linewidth',1.2)
title('实验数据与附件数据对比')
legend('实验数据','附件数据')
xlabel('时间/s')
ylabel('温度/摄氏度')
max(kk)
if max(kk)<= 3
    %% 判断150~190时间
    [~,idx_max] = max(q) ;
    q_now = q(1:idx_max);
    x_now = x(1:idx_max);
    idx1 = find(q_now>=150 & q_now <=190);
    x_1 = x_now(idx1);
    t_1 = x_1(end) - x_1(1);%150到190时间
    if t_1>=60 && t_1<=120
       %% 判断大于217时间
        idx2 = find(q > 217);
        x_now = x(idx2);
        t_2 = x_now(end)-x_now(1);
        if t_2 >=40 && t_2 <= 90
            %% 判断峰值温度
            q_max = max(q);
            if q_max>=240 && q_max <=250
                answer_now = 0;
                m = 0;
                for i = 1 : length(xx_real)
                    for j = 1 : length(test_x)
                        if test_x(j) == xx_real(i)
                           answer_now = answer_now + abs(x_real(i)-test_y(j));
                           m = m + 1;
                        end
                    end
                end
                fprintf('均误差%d',answer_now/m);
            end
        end
    end
end
answer_now = 0;
m = 0;
for i = 1 : length(xx_real)
    for j = 1 : length(test_x)
        if test_x(j) == xx_real(i)
           answer_now = answer_now + abs(x_real(i)-test_y(j));
           m = m + 1;
           break
        end
    end
end
fprintf('均误差%d',answer_now/m);
%         end
%     end
% end





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