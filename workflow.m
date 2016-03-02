% workflow

m = 1; b = 1; k = 1; xd = 2;
A = @(t)[0 1 ; -k/m -b/m];
B = @(t)[0; 1/m];
Q = @(t)eye(2);
R = @(t)1;
S = eye(2);
[K,P] = dare(A,B,Q,R,S,5,100);
f = @(t,x,u)A(t)*x+B(t)*u(1);
u = @(x,t)-K(t)*x-inv([1 0]*inv(A(t)-B(t)*K(t))*B(t))*xd;
[X,U,t] = dynSim(f,u,[-1 ; 0],5,0.1);
subplot(2,1,1); plot(t,X); xlabel('t');ylabel('x');
subplot(2,1,2); plot(t,U); xlabel('t');ylabel('u');