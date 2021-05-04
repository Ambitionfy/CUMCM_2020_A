function [nnn,qqq,fitness] = fun2(u_1,u_2,u_3,u_4,v,ke1,ks)
global rol T h ms lemda
fitness = 100000000;
sk = 0;
answer = [];%样本点
ke = 144;
% ks = 9.3;
[p,p8,p9] = get_kall(ke,ks,u_1,u_2,u_3,u_4);
% [p] = get_keks(ke,ks);
% p8 = get_p8(ke,ks);
% p9 = get_p9(ke,ks);
f1 = @(x)p(1,1)*x^4+p(1,2)*x^3+p(1,3)*x^2+p(1,4)*x+p(1,5);%175-175
f2 = @(x)p(2,1)*x^4+p(2,2)*x^3+p(2,3)*x^2+p(2,4)*x+p(2,5);%175-195
f3 = @(x)p(3,1)*x^4+p(3,2)*x^3+p(3,3)*x^2+p(3,4)*x+p(3,5);%195-235
f4 = @(x)p(4,1)*x^4+p(4,2)*x^3+p(4,3)*x^2+p(4,4)*x+p(4,5);%235-255
f5 = @(x)p(5,1)*x^4+p(5,2)*x^3+p(5,3)*x^2+p(5,4)*x+p(5,5);%255-255
% 炉前温度
f8 = @(x)p8(1)*x^9+p8(2)*x^8+p8(3)*x^7+p8(4)*x^6+p8(5)*x^5+p8(6)*x^4+p8(7)*x^3+p8(8)*x^2+p8(9)*x+p8(10);
f9 = @(x)p9(1)*x^9+p9(2)*x^8+p9(3)*x^7+p9(4)*x^6+p9(5)*x^5+p9(6)*x^4+p9(7)*x^3+p9(8)*x^2+p9(9)*x+p9(10);

L = 50+30.5*11+5*10;
% v = 78/60;%cm/s
t_real = L/v;
h = 1e-4;
% ke1 = 32.68;
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
        u0 = u_1;
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(2) && i<=mm_real(3))
        %空隙1
        l_now = (i-mm_real(2))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(3) && i<=mm_real(4))
        %温区2
        u0 = u_1;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(4) && i<=mm_real(5))
        %空隙2
        l_now = (i-mm_real(4))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(5) && i<=mm_real(6))
        %温区3
        u0 = u_1;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(6) && i<=mm_real(7))
        %空隙3
        l_now = (i-mm_real(6))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(7) && i<=mm_real(8))
        %温区4
        u0 = u_1;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(8) && i<=mm_real(9))
        %空隙4
        l_now = (i-mm_real(8))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(9) && i<=mm_real(10))
        %温区5
        u0 = u_1;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(10) && i<=mm_real(11))
        %空隙5
        l_now = (i-mm_real(10))/100*v;
        u0 = f2(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(11) && i<=mm_real(12))
        %温区6
        u0 = u_2;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(12) && i<=mm_real(13))
        %空隙6
        l_now = (i-mm_real(12))/100*v;
        u0 = f3(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(13) && i<=mm_real(14))
        %温区7
        u0 = u_3;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(14) && i<=mm_real(15))
        %空隙7
        l_now = (i-mm_real(14))/100*v;
        u0 = f4(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(15) && i<=mm_real(16))
        %温区8
        u0 = u_4;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(16) && i<=mm_real(17))
        %空隙8
        l_now = (i-mm_real(16))/100*v;
        u0 = f5(l_now);
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
    elseif (i>mm_real(17) && i<=mm_real(18))
        %温区9
        u0 = u_4;
        A(1,1:2) = [h*ke1/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke1/lemda+1];
        f(1) = h*ke1*u0/lemda;
        f(ms+1) = h*ke1*u0/lemda;
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
        l_now = (i-mm_real(18))/100*v;
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
x = x/v;
x(idx) = [];
nnn = x;
qqq = q;
%% 1
% plot(x,q,'b-','linewidth',1.2,'k-.')
% hold on
% plot(x,217*ones(size(x)),'k-.','linewidth',1.2)
% title('当前条件下炉温曲线')

%% 判断
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
for i = 2 : length(x_real)
    kk = [kk abs(x_real(i)-x_real(i-1))/abs(xx_real(i)-xx_real(i-1))];
end
if max(kk)<= 3
    %% 判断150~190时间
    [~,idx_max] = max(q) ;
    q_now = q(1:idx_max);
    x_now = x(1:idx_max);
    idx1 = find(q_now>=150 & q_now <=190);
    x_1 = x_now(idx1);
    t_1 = x_1(end) - x_1(1)%150到190时间
               fitness = 0;
               idx = find(q>217);
               [~,idx1] = max(q);
               q_real = q(idx:idx1);
               x_real = x(idx:idx1);
               for i = 1 : length(x_real)-1
                   delta_t = abs(x_real(i+1)-x_real(i));
                   fitness = fitness + (q_real(i)-217 + q_real(i+1)-217)/2*delta_t;
               end
               if fitness <= 8000
               idx = find(q>217);
               q_real = q(idx);
               x_real = x(idx);
               y = floor(q_real-217);
               for i = 1 : length(x_real)
                   answer = [answer x_real(i)*ones(1,y(i))];
               end
               %计算样本点的SK
               answer_mean = mean(answer);
               sigma = std(answer);
               n = length(answer);
               for i = 1 : length(answer)
                  sk = sk + (answer(i)-answer_mean)^3; 
               end
               sk = sk/(sigma^3*n);
               fitness = sk;
                else
                fitness = 1e6;
              
               end
end


function y = get_c(T)
%传入当前温度，获得当前芯片比热容
%比热容-温度多项式系数
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