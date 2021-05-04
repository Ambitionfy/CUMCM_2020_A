function x=chase(a,b,c,f) 
%求解线性方程组 Ax=f, 其中 A是三对角阵
%a是矩阵 A的下对角线元素 a(1)=0
%b是矩阵 A的对角线元素
%c是矩阵 A的上对角线元素 c(n)=0
%f是方程组的右端向量
n=length(f); 
x=zeros(1,n);
y=zeros(1,n); 
d=zeros(1,n);
z= zeros(1,n); 
%预处理
d(1)=b(1); 
for i=1:n-1 
z(i)=c(i)/d(i); 
d(i+1)=b(i+1)-a(i+1)*z(i); 
end
%追的过程
y(1)=f(1)/d(1); 
for i=2:n 
 y(i)=(f(i)-a(i)*y(i-1))/d(i); 
end
%赶的过程
x(n)=y(n); 
for i=n-1:-1:1 
x(i)=y(i)-z(i)*x(i+1); 
end