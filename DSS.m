function [x_traj,u_traj,t,J] = DSS(L,M,h,r,f,x_0,m,T,N)

% initialize
w_0 = ones(m*(N+1),1)*0.1;

% defining objective function
fun = @(w)bolzaDSS(w,L,M,f,x_0,T,N);

% defining constraints
nonlcon = @(w)constraintDSS(w,f,h,r,x_0,T,N);

% solving transformed NLP
options = optimoptions('fmincon','MaxFunEvals',20000,'UseParallel',1,'Algorithm','sqp', 'TolCon',1.0000e-01,'TolFun', 1.0000e-01,'TolX',0.1);
[w,J] = fmincon(fun,w_0,[],[],[],[],[],[],nonlcon,options);
U = reshape(w,numel(w)/(N+1),N+1);
u_traj = @(t)getControlTraj(t,U,T);
% fill in this part
[X,t] = forSim(f,x_0,@(x,t) u_traj(t),T,N);
x_traj = @(t)getStateTraj(t,X,T);

end

function J = bolzaDSS(w,L,M,f,x_0,T,N)

U_guess = reshape(w,numel(w)/(N+1),N+1);
u_guess_traj = @(t)getControlTraj(t,U_guess,T);

% forward simulate for X
% fill in this part
[X_sim,~] = forSim(f, x_0, @(x,t) u_guess_traj(t), T, N);
x_sim_traj = @(t)getStateTraj(t,X_sim,T);
% numerically integrate the Lagrange term
n = numel(x_0);
u_a = @(Lambda,t)augmentControl(x_sim_traj,u_guess_traj,t);
L_a = @(Lambda,u_a,t)augmentLagrange(L,u_a,t,n);

Lambda_0 = 0;

Lambda = forSim(L_a,Lambda_0,u_a,T,N);

% evaluate total Bolza objective cost
% fill in this part
J = Lambda(end) + M(X_sim(end),T);

% augmented Lagrange term to work with forSim
    function L_a = augmentLagrange(L,u_a,t,n)
        x = u_a(1:n);
        u = u_a(n+1:end);
        L_a = L(x,u,t);
    end

    function u_a = augmentControl(x,u,t)
        u_a = [x(t);u(t)];
    end

end

function [inCon,eqCon] = constraintDSS(w,f,h,r,x_0,T,N)
% constraints for fmincon
U_guess = reshape(w,numel(w)/(N+1),N+1);
u_guess_traj = @(t)getControlTraj(t,U_guess,T);

% fill in this part
[X_sim,~] = forSim(f, x_0, @(x,t) u_guess_traj(t), T, N);

% allocate space for inCOn
n_h = size(h(X_sim(:,1),U_guess(:,1)),1); % total number of discretized path constraints
inCon = zeros(n_h*(N+1),1);

% fill in the inequality constrants
for i = 1:N+1,
    % fill in this part
    inCon(n_h*(i-1)+1:n_h*i) = h(X_sim(:,i), U_guess(:,i));
end

% fill in this part
inCon(end+1:end+size(r(X_sim(:,end),T),1)) = r(X_sim(:,end),T);

% terminal constraints
% fill in this part
eqCon = [];
%eqCon(end+1:end+size(r(X_sim(:,end),T),1)) = r(X_sim(:,end),T);

end

function x = getStateTraj(t,X,T)
%   returns a state value evaluated at time t based on
%   parameterized value X and final time T.
%   This approach uses piece-wise linear
n = size(X,1);
N = size(X,2)-1;
numT = length(t);
x = zeros(n,numT);
for i = 1:length(t)
    if t(i)<0
        x(:,i) = X(:,1);
    elseif t(i)>=T
        x(:,i) = X(:,end);
    else
        time_vector = linspace(0,T,N);
        idx =  floor(t(i)*N/T)+1;
        dt = T/N;
        x(:,i) = (X(:,idx+1)-X(:,idx))*(t(i)-time_vector(idx))/dt+X(:,idx);
    end
end
end

function u = getControlTraj(t,U,T)
%   returns a control input value evaluated at time t based on
%   parameterized value U and final time T.
%   This approach uses piece-wise constant

m = size(U,1);
N = size(U,2)-1;
numT = length(t);
u = zeros(m,numT);
for i = 1:length(t)
    if t(i)<0
        u(:,i) = U(:,1);
    elseif t(i)>T
        u(:,i) = U(:,end);
    else
        idx =  floor(t(i)*N/T)+1;
        u(:,i) = U(:,idx);
    end
end

end