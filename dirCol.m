function [x_traj,u_traj,t,J] = dirCol(L,M,h,r,f,x_0,m,T,N)
%%DIRCOL - generates a set of optimal trajectories based on given CTOCP
%%formulation with direct collocation method
%   [X,U,T,J] = DIRCOL(L,M,H,R,F,X0,MU,TF,N) generates optimal state
%   trajectory X and control input trajectory U that minimize the bolza
%   cost J of the following form
%                J = int(L(X,U,t),t=0->TF)+M(X(TF),TF)
%   The trajectories are subjected to path constraint H and terminal
%   constraint R. F is the dynamics of the states. X0 is the initial
%   states. MU is the number of control input. TF is final time. N is the
%   number of discretized interval.
%   
%   X and U are function handles in term of time
%   L is a function handle in term of numeric states, numeric control
%   input, and time.
%   M is a function handle in term of numeric final state, and final time.
%   R is a function handle in term of numeric final state, and final time.
%   F is a function handle in term of numeric states, numeric control
%   input, and time.
%   **** the terminal constraint is an equality *****
% Author: Tucker Haydon
% Modified by Pi Thanacha Choopojcharoen

% initialize
n = numel(x_0);
w_0 = 0.00001*ones((n+m)*(N+1) + N,1);

% defining objective function
fun = @(w)bolzaDircol(w,L,M,x_0,m,T,N);

% defining constraints
nonlcon = @(w)conDircol(w,x_0,m,T,N,f,r,h);

% solving transformed NLP
options = optimoptions('fmincon','MaxFunEvals',1000000,'MaxIter',100000,'UseParallel',1);
[w,J] = fmincon(fun,w_0,[],[],[],[],[],[],nonlcon,options);

% Define the time vector
t = linspace(0,T,N+1);

% extract U and X from w
X = reshape(w(1:n*(N+1)),n,N+1);
U = reshape(w(n*(N+1)+1:(n+m)*(N+1)),m,N+1);
U = [U(:,1:end-1) U(:,end-1)];

% Convert discretized trajectories to continuous state and control trajectories
u_traj = @(t)getControlTraj(t,U,T);
x_traj = @(t)getStateTraj(t,X,U,T,f);

end
function J = bolzaDircol(w,L,M,x_0,m,T,N)

% Constants
n = numel(x_0);
t = linspace(0,T,N+1);

% Extract X and U
X = reshape(w(1:n*(N+1)),n,N+1);
U = reshape(w(n*(N+1)+1:(n+m)*(N+1)),m,N+1);
U = [U(:,1:end-1) U(:,end-1)];

% The final cost
J = sum(L(X, U, t)) + M(X(:,end),T);


end
function [inCon, eqCon] = conDircol(w,x_0,m,T,N,f,r,h)

% Constants
n = numel(x_0);
t = linspace(0,T,N+1);
dt = T/N;
t_c = t(1:end-1) + dt/2;

% Extract X and U
X = reshape(w(1:n*(N+1)),n,N+1);
U = reshape(w(n*(N+1)+1:(n+m)*(N+1)),m,N+1);
U = [U(:,1:end-1) U(:,end-1)];

% Define B
B = [zeros(n) 1/dt*eye(n) 1/(dt.^2)*eye(n) 3/(4*(dt.^3))*eye(n)];

% Define A
A = [eye(n) 1/2*eye(n) 1/4*eye(n) 1/8*eye(n)];

% Define C
C  = [eye(n) zeros(n) zeros(n) zeros(n);
    eye(n) eye(n) eye(n) eye(n);
    zeros(n) 1/dt*eye(n) zeros(n) zeros(n);
    zeros(n) 1/dt*eye(n) 2/(dt.^2)*eye(n) 3/(dt.^3)*eye(n)];

% Define D
D = [X(:,1:end-1);
    X(:,2:end);
    f(X(:,1:end-1), U(:,1:end-1), t(1:end-1));
    f(X(:,2:end),U(:,2:end), t(2:end))];

% Calculate P and dP: The collocation point and derivative
V = C\D;
P = A*V; % Value of the polynomial at collocation points
dP = B*V; % Time derivative of polynimial at collocation points


row = size(X,1);
col = size(X,2);
U_c = getControlTraj(t_c,U,T);

% Equality constraints
% Boundary constraint!
% Colloction constraint!
% Terminal contraint! 
eqCon = [X(:,1) - x_0;
    reshape(f(P,U_c,t),row*N,1) - reshape(dP,row*N,1)];

% Inequality constraints
% Path constraint

inCon = [reshape(h(X, U),size(h(X,U),1)*size(h(X,U),2),1);
    r(X(:,end),T)];
% inCon = [];

end
function u = getControlTraj(t,U,T)
%   returns a state value evaluated at time t based on
%   parameterized value X and final time T.
%   This approach uses piece-wise linear
n = size(U,1);
N = size(U,2)-1;
numT = length(t);
u = zeros(n,numT);
for i = 1:length(t)
    if t(i)<0
        u(:,i) = U(:,1);
    elseif t(i)>=T
        u(:,i) = U(:,end);
    else
        time_vector = linspace(0,T,N);
        idu =  floor(t(i)*N/T)+1;
        dt = T/N;
        u(:,i) = (U(:,idu+1)-U(:,idu))*(t(i)-time_vector(idu))/dt+U(:,idu);
    end
end
end
function x = getStateTraj(t,X,U,T,f)
%   returns a state value evaluated at time t based on
%   parameterized value X and final time T.
n = size(X,1);
N = size(X,2)-1;
numT = length(t);
x = zeros(n,numT);

dt = T/N;
U = [U(:,1:end-1) U(:,end-1)];

% Define C
C  = [eye(n) zeros(n) zeros(n) zeros(n);
    eye(n) eye(n) eye(n) eye(n);
    zeros(n) 1/dt*eye(n) zeros(n) zeros(n);
    zeros(n) 1/dt*eye(n) 2/(dt.^2)*eye(n) 3/(dt.^3)*eye(n)];

% Define D
D = [X(:,1:end-1);
    X(:,2:end);
    f(X(:,1:end-1), U(:,1:end-1), t(1:end-1));
    f(X(:,2:end),U(:,2:end), t(2:end))];

V = C\D;

for i = 1:length(t)
    if t(i)<0
        x(:,i) = X(:,1);
    elseif t(i)>=T
        x(:,i) = X(:,end);
    else
        time_vector = linspace(0,T,N);
        idx =  floor(t(i)*N/T)+1;
        dt = T/N;
        s = (t(i)-time_vector(idx))/dt;
        A = [eye(n) s*eye(n) s^2*eye(n) s^3*eye(n)];
        x(:,i) = A*V(:,idx);
    end
end

end