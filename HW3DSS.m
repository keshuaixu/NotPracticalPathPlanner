clc;

L = @(x,u,t)(u^2);
M = @(x,T)(0);
u_max = 10;
h = @(x,u)[u-u_max;-u-u_max];
e = 0.0001; % error 
r = @(x,T) [-e-((1-cos(x(1)))/2);((1-cos(x(1)))/2)-e];


m = 0.1;
l = 1;
g = 9.81;
b = 0.1;

f = @(x,u,t)[x(2);(m*g*l*sin(x(1))-b*x(2)+u)/(m*l^2)];
x_0 = [ pi/2 ; 0 ];
T = 1;
N = 30;

m = 1; % number of control input
tic;
[x,u,t,J] = DSS(L,M,h,r,f,x_0,m,T,N);
toc
%%

X = x(t);
U = u(t);
subplot(2,1,1)
plot(t,X)
subplot(2,1,2)
plot(t,U)

figure
hold on;
t_previous = -T/N;
for i = 1:size(X,2)
    
    p = plot([0 l*sin(X(1,i))],[0 l*cos(X(1,i))],'b','LineWidth',3);
    axis([-1 1 -1 1])
    axis equal
    
    pause(t(i)-t_previous);
    t_previous = t(i);
    if i ~= size(X,2)
        delete(p);
        
    end
end