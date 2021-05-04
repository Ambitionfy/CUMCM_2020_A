function test8
%�ڶ���
clc,clear,close all
global rol T h ms lemda
ke = 144;
ks = 10.5;
[p,p8,p9] = get_kall(ke,ks,182,203,237,254);
% [p] = get_keks(ke,ks);
% p6 = get_p6(ke,ks);
% p8 = get_p8(ke,ks);
% p9 = get_p9(ke,ks);
f1 = @(x)p(1,1).*x.^4+p(1,2)*x.^3+p(1,3)*x.^2+p(1,4).*x+p(1,5);%175-175
f2 = @(x)p(2,1).*x.^4+p(2,2)*x.^3+p(2,3)*x.^2+p(2,4).*x+p(2,5);%175-195
f3 = @(x)p(3,1).*x.^4+p(3,2)*x.^3+p(3,3)*x.^2+p(3,4).*x+p(3,5);%195-235
f4 = @(x)p(4,1).*x.^4+p(4,2)*x.^3+p(4,3)*x.^2+p(4,4).*x+p(4,5);%235-255
f5 = @(x)p(5,1).*x.^4+p(5,2)*x.^3+p(5,3)*x.^2+p(5,4).*x+p(5,5);%255-255
f6 = @(x)p6(1)*x.^9+p6(2)*x.^8+p6(3)*x.^7+p6(4)*x.^6+p6(5)*x.^5+p6(6)*x.^4+p6(7)*x.^3+p6(8)*x.^2+p6(9)*x+p6(10);
% ¯ǰ�¶�
f8 = @(x)p8(1)*x.^9+p8(2)*x.^8+p8(3)*x.^7+p8(4)*x.^6+p8(5)*x.^5+p8(6)*x.^4+p8(7)*x.^3+p8(8)*x.^2+p8(9)*x+p8(10);
f9 = @(x)p9(1)*x.^9+p9(2)*x.^8+p9(3)*x.^7+p9(4)*x.^6+p9(5)*x.^5+p9(6)*x.^4+p9(7)*x.^3+p9(8)*x.^2+p9(9)*x+p9(10);
L = 50+30.5*11+5*10;
result1 = 100:-1:65;
result2 = [];
for v = 100/60 : -1/60 : 65/60
t_real = L/v;
h = 1e-4;
% ke = 29.67;
ke = 32.68;
T = 0.01;
Tmax = fix(t_real/T);
lemda = 98.9;
rol = 2330;
c = get_c(25+273);%��ʼ������
k = lemda/(rol*c);
m=0.15;%ʵ�ʳ���
ms = m*100;
mm(1) = 25;
for i = 1 : 10
    mm = [mm 30.5 5];
end
mm = [mm 30.5 25];%�ܾ���
mm_fix = mm/v*100;
mm_real(1)=mm_fix(1);
for i = 2 : length(mm_fix)
    mm_real(i) = mm_fix(i)+mm_real(i-1);
end
mm_real = fix(mm_real);%��ÿ��ʱ��
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


a = zeros(1,ms+1); %�¶Խ�
b = zeros(1,ms+1); %�жԽ�
c = zeros(1,ms+1); %�϶Խ�
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
        %¯ǰ
        l_now = i/100*v;
        u0 = f8(l_now);
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ks*u0/lemda;
    elseif (i>mm_real(1) && i<=mm_real(2))
        %����1
        u0 = 175;
        A(1,1:2) = [h*ks/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ks*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(2) && i<=mm_real(3))
        %��϶1
        l_now = (i-mm_real(2))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(3) && i<=mm_real(4))
        %����2
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(4) && i<=mm_real(5))
        %��϶2
        l_now = (i-mm_real(4))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(5) && i<=mm_real(6))
        %����3
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(6) && i<=mm_real(7))
        %��϶3
        l_now = (i-mm_real(6))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(7) && i<=mm_real(8))
        %����4
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(8) && i<=mm_real(9))
        %��϶4
        l_now = (i-mm_real(8))/100*v;
        u0 = f1(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(9) && i<=mm_real(10))
        %����5
        u0 = 175;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(10) && i<=mm_real(11))
        %��϶5
        l_now = (i-mm_real(10))/100*v;
        u0 = f2(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(11) && i<=mm_real(12))
        %����6
        u0 = 195;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(12) && i<=mm_real(13))
        %��϶6
        l_now = (i-mm_real(12))/100*v;
        u0 = f3(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(13) && i<=mm_real(14))
        %����7
        u0 = 235;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(14) && i<=mm_real(15))
        %��϶7
        l_now = (i-mm_real(14))/100*v;
        u0 = f4(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(15) && i<=mm_real(16))
        %����8
        u0 = 255;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(16) && i<=mm_real(17))
        %��϶8
        l_now = (i-mm_real(16))/100*v;
        u0 = f5(l_now);
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
    elseif (i>mm_real(17) && i<=mm_real(18))
        %����9
        u0 = 255;
        A(1,1:2) = [h*ke/lemda+1 -1];
        A(ms+1,ms:ms+1) = [ -1 h*ke/lemda+1];
        f(1) = h*ke*u0/lemda;
        f(ms+1) = h*ke*u0/lemda;
%     elseif (i>mm_real(18) && i<=mm_real(19))
%         %��϶9
%         l_now = (i-mm_real(18))/100*v;
%         u0 = f6(l_now);
%         A(1,1:2) = [h*ke/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ke*255/lemda;
%         f(ms+1) = h*ks*u0/lemda;
%     elseif (i>mm_real(19) && i<=mm_real(20))
%         %����10
%         u1 = 25;
% %         l_now = (i-mm_real(19))/100*70/60;
% %         u0 = f7(l_now);
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u0/lemda;
%         f(ms+1) = h*ks*u1/lemda;
%     elseif (i>mm_real(20) && i<=mm_real(21))
%         %��϶10
%         u0 = 25;
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u0/lemda;
%         f(ms+1) = h*ks*u0/lemda;
%     elseif (i>mm_real(21) && i<=mm_real(22))
%         %����11
%         u0 = 25;
%         A(1,1:2) = [h*ks/lemda+1 -1];
%         A(ms+1,ms:ms+1) = [ -1 h*ks/lemda+1];
%         f(1) = h*ks*u0/lemda;
%         f(ms+1) = h*ks*u0/lemda;
%     elseif (i>mm_real(22))
%         %¯��
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
    c_now = get_c(u0+273);%��ǰʱ�̱�����
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
% title('�¶���ʱ�䡢�ռ�ֲ�ͼ')
% xlabel('�ռ�/0.1mm')
% ylabel('ʱ��/s')
% zlabel('�¶�/���϶�')
% figure(2)
x = v/100:v/100:435.5+v/100;
x(idx) = [];
x = x/v;
% axis([0 450 0 250])
%% �ж�б��
xmin = ceil(x(1));
xmax = floor(x(end));
idx = [];
for i = xmin:xmax
    idx = [idx find(x==i)];
end
x_real = q(idx);
xx_real = x(idx);
kk = [];%б��
for i = 2 : min(length(x),length(q))
    kk = [kk abs(q(i)-q(i-1))/abs(x(i)-x(i-1))];
end
% [~,idx] = max(kk);
% hold on
% plot(x(idx),q(idx),'*')
% plot(x,q)
if max(kk)<= 3
    %% �ж�150~190ʱ��
    [~,idx_max] = max(q) ;
    q_now = q(1:idx_max);
    x_now = x(1:idx_max);
    idx1 = find(q_now>=150 & q_now <=190);
    x_1 = x_now(idx1);
    t_1 = x_1(end) - x_1(1);%150��190ʱ��
    if t_1>=60 && t_1<=120
       %% �жϴ���217ʱ��
        idx2 = find(q > 217);
        x_now = x(idx2);
        t_2 = x_now(end)-x_now(1);
        if t_2 >=40 && t_2 <= 90
            %% �жϷ�ֵ�¶�
            q_max = max(q);
            if q_max>=240 && q_max <=250
                t = 1;
                result2 = [result2 1];
            else
                result2 = [result2 0];
            end
        else
            result2 = [result2 0];
        end
    else
        result2 = [result2 0];
    end
else
    result2 = [result2 0];
end
end
idx = find(result2 == 1);
result1(idx)
figure
plot(result1,result2,'*','linewidth',4)
title('��ͬ�ٶ����������')
xlabel('�ٶ�/cm/min')
ylabel('�������')




function y = get_c(T)
%���뵱ǰ�¶ȣ���õ�ǰоƬ������
%������-�¶ȶ���ʽϵ��
% x = [100 200 300 400 600 800 1000 1200 1500];
% y = [259 556 712 790 867 913 946 967 992];
% p = polyfit(x,y,6);
% y = polyval(p,T);
p_t = 1e2*[-0.000000000000000   0.000000000000111  -0.000000000252512   0.000000299838736  -0.000197353916376   0.071405959681669  -2.853032199050470];
y = polyval(p_t,T);

function [a,b,c] = get_abc(cc,A)
global rol T h ms lemda
%��������ݣ����abc
k = lemda/(rol*cc);
r = T * k /  h^2;
for i = 2 : ms
   A(i,i-1) = -r;
   A(i,i) = 1+2*r;
   A(i,i+1) = -r;
end
a = zeros(1,ms+1); %�¶Խ�
b = zeros(1,ms+1); %�жԽ�
c = zeros(1,ms+1); %�϶Խ�
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

    
    