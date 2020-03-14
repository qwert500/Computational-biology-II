%-------------------------
%      Comp. Bio II
%     Example sheet 4
%         Task 2a
%   Deterministic approx.        
%-------------------------
clc, clear
%-------------------------
%     Initialization
%-------------------------
m=0.01;
s=0.1;
x1=0;
x2=0.0025;
limit=10^(-6);
epsilon1=1;
epsilon2=1;

w1AA=1-s;
w1Aa=1-s/2;
w1aa=1;
w2AA=1;
w2Aa=1-s/2;
w2aa=1-s;
t=0;
while epsilon1>limit || epsilon2>limit
  t=t+1;
  w1(t)=(x1(t)^2*w1AA+2*x1(t)*(1-x1(t))*w1Aa+(1-x1(t))^2*w1aa)*(1-m)+...
    ((x2(t))^2*w1AA+2*x2(t)*(1-x2(t))*w1Aa+(1-x2(t))^2*w1aa)*m;
  
  w2(t)=(x2(t)^2*w2AA+2*x2(t)*(1-x2(t))*w2Aa+(1-x2(t))^2*w2aa)*(1-m)+...
    ((x1(t))^2*w2AA+2*x1(t)*(1-x1(t))*w2Aa+(1-x1(t))^2*w2aa)*m;
  
  x1(t+1)=((x1(t))^2*w1AA+x1(t)*(1-x1(t))*w1Aa)*(1-m)/(w1(t))+...
    ((x2(t))^2*w1AA+x2(t)*(1-x2(t))*w1Aa)*m/(w1(t));
  
  x2(t+1)=((x2(t))^2*w2AA+x2(t)*(1-x2(t))*w2Aa)*(1-m)/(w2(t))+...
    ((x1(t))^2*w2AA+x1(t)*(1-x1(t))*w2Aa)*m/(w2(t));
  
  epsilon1=x1(t+1)-x1(t);
  epsilon2=x2(t+1)-x2(t);
end
hold on
plot(x1,x2)
plot(x1(end),x2(end),'o')
x1(end)
x2(end)
%%
theory_x1m01s01=x1(end);
theory_x2m01s01=x2(end);












