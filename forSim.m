function [X,t] = forSim(f,x_0,u,T_f,N)
%FORSIM - Forward Simulate the dynamics system
%   [X,T] = FORSIM(F,X0,U,TF,N) forward simulates a given dynamic system F
%
%   x'(t) = f(x(t),u(x(t),t),t)
%
%   with a given initial states X0, starting at time t = 0 and ending at 
%   t = TF. U is a control policy which is a function handle in term of
%   state variable and time.

t = linspace(0,T_f,N+1);
X = zeros(numel(x_0),N+1);
dt = T_f/N;
X(:,1) = x_0;
for i = 1:N,
    X(:,i+1) = rk4(f,X(:,i),u,t(i),dt);
end

end

function X_next = rk4(f,X_i,u,t_i,dt)
% rk4
k1 = f(X_i,u(X_i,t_i),t_i);
k2 = f(X_i+k1*dt/2,u(X_i+k1*dt/2,t_i+dt/2),t_i+dt/2);
k3 = f(X_i+k2*dt/2,u(X_i+k2*dt/2,t_i+dt/2),t_i+dt/2);
k4 = f(X_i+k3*dt,u(X_i+k3*dt,t_i+dt),t_i+dt);
X_next = X_i + (k1+2*k2+2*k3+k4)*dt/6;

end
