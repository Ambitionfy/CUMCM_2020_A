function x=chase(a,b,c,f) 
%������Է����� Ax=f, ���� A�����Խ���
%a�Ǿ��� A���¶Խ���Ԫ�� a(1)=0
%b�Ǿ��� A�ĶԽ���Ԫ��
%c�Ǿ��� A���϶Խ���Ԫ�� c(n)=0
%f�Ƿ�������Ҷ�����
n=length(f); 
x=zeros(1,n);
y=zeros(1,n); 
d=zeros(1,n);
z= zeros(1,n); 
%Ԥ����
d(1)=b(1); 
for i=1:n-1 
z(i)=c(i)/d(i); 
d(i+1)=b(i+1)-a(i+1)*z(i); 
end
%׷�Ĺ���
y(1)=f(1)/d(1); 
for i=2:n 
 y(i)=(f(i)-a(i)*y(i-1))/d(i); 
end
%�ϵĹ���
x(n)=y(n); 
for i=n-1:-1:1 
x(i)=y(i)-z(i)*x(i+1); 
end