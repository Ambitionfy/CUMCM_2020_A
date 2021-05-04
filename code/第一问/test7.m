function test7
%功能——输出第一问间隔温度
ke = 144;
ks = 10.5;
[p,p8,p9] = get_kall(ke,ks,173,198,230,257);
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
% 炉前温度
f8 = @(x)p8(1)*x.^9+p8(2)*x.^8+p8(3)*x.^7+p8(4)*x.^6+p8(5)*x.^5+p8(6)*x.^4+p8(7)*x.^3+p8(8)*x.^2+p8(9)*x+p8(10);
f9 = @(x)p9(1)*x.^9+p9(2)*x.^8+p9(3)*x.^7+p9(4)*x.^6+p9(5)*x.^5+p9(6)*x.^4+p9(7)*x.^3+p9(8)*x.^2+p9(9)*x+p9(10);

%% 计算
result = [];
x1 = 0.01 : 0.01 : 25;
x2  =0.01 : 0.01 : 5;
x3  = 0.01 : 0.01:30.5;
x4 = 0.01 : 0.01  : 40.5;
result = [f8(x1) 173*ones(size(x3))  f1(x2) 173*ones(size(x3))  f1(x2) 173*ones(size(x3))  f1(x2)  173*ones(size(x3))  f1(x2) 173*ones(size(x3))  f2(x2)  198*ones(size(x3)) f3(x2) ];
% result = [result 230*ones(size(x3)) f4(x2) 255*ones(size(x3))  f5(x2)  255*ones(size(x3)) f6(x2)  25*ones(size(x3))  25*ones(size(x2))  25*ones(size(x3))  25*ones(size(x1))];
result = [result 230*ones(size(x3)) f4(x2) 257*ones(size(x3))  f5(x2)  257*ones(size(x3)) f9(x4)  25*ones(size(x3))  25*ones(size(x1))];
x = 0.01 : 0.01:435.5;
length(x)
length(result)
plot(x,result,'k','linewidth',1.2)
title('回焊炉内部温度实际图');
xlabel('空间/cm')
ylabel('温度/摄氏度')