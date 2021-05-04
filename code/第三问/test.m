function test
x1 = 1 : 2;
x2 = 3:30;
y = [561*ones(size(x1)),456.6368*ones(size(x2))];
x = [x1 x2];
plot(x,y,'k','linewidth',1.2)
title('最优个体适应度')
xlabel('迭代次数')
ylabel('适应度')
% x1 = 1 : 11;
% x2 = 12:30;
% y = [0.243*ones(size(x1)),0.1419*ones(size(x2))];
% x = [x1 x2];
% plot(x,y,'k','linewidth',1.2)
% title('最优个体适应度')
% xlabel('迭代次数')
% ylabel('适应度')