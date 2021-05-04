function test
x = [178.1777,184.5079,235,261,87/60,22,10.5;178.1777,184.5079,235,261,87/60,29,10.5;178.1777,184.5079,235,261,87/60,25,3.8];
figure(1)
hold on
[x1,y1,z1] = fun2(x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),x(1,6),x(1,7)) ;
[x2,y2,z2] = fun2(x(2,1),x(2,2),x(2,3),x(2,4),x(2,5),x(2,6),x(2,7)) ;
[x3,y3,z3] = fun2(x(3,1),x(3,2),x(3,3),x(3,4),x(3,5),x(3,6),x(3,7)) ;
plot(x1,y1,'k--','linewidth',1.2)
plot(x2,y2,'k-','linewidth',1.2)
plot(x3,y3,'k-.','linewidth',1.2)
x4 = x1(1) : 10 : x1(end);
plot(x4,217*ones(size(x4)),'k-o','linewidth',1.2)
title('不同SK下的炉温曲线局部图')
legend('SK = -0.0588','SK = -0.1631','SK = 0.0703','217温度线')