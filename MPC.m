function u = MPC(x_k,t_k,x_r,u_r,f,h,Q,R,T_h,N)
% u is a function handle in term of time t
R_k = R(t_k);
Q_k = Q(t_k);
S_k = S(t_k);
w_0 = zeros(size(R_k),N+1);
U = fmincon(@(w)objFun(w,x_k,t_k,x_r,u_r,f,Q_k,R_k,S_k,T_h,N),w_0,[],[],[],[],[],[],@(u)cons(w,x_k,t_k,f,h,R_k,T_h,N));
u_optimal = reshape(U,size(R_k,1),N+1);
u = @(t)control(t,t_k,u_optimal,T_h,N);
u = @(t)u(t);
    function ut = control(t,t_k,u_optimal,T_h,N)
        i = floor(t*N/T_h)+1;
        ut = u_optimal(:,i)+(u_optimal(:,i+1)-u_optimal(:,i))*(t-t_k)/T_h;
    end
end

function J = objFun(w,x_k,t_k,x_r,u_r,f,Q_k,R_k,S_k,T_h,N)
u_guess = reshape(w,size(R_k,1),N+1);
u = @(x,t)control(t,t_k,u_guess,T_h,N);
    function ut = control(t,t_k,u_guess,T_h,N)
        i = floor(t*N/T_h)+1;
        ut = u_guess(:,i)+(u_guess(:,i+1)-u_guess(:,i))*(t-t_k)/T_h;
    end
x_sim = dynSim(f,u,x_k,T,dt);
t = (0:T_h/N:T_h)+t_k;
e_x = e_sim-x_r(t); % discretized state error trajectory
e_u = u_guess-u_r(t); % discretized control error trajectory

e_X = reshape(e_x(:,1:end-1),size(Q_k,1)*N,1);
e_U = reshape(e_u(:,1:end-1),size(R_k,1)*N,1);

intL = e_X'*Q_k*e_X + e_U'*R_k*e_U;
M = x_sim(:,end)'*S_k*x_sim(:,end);
J = intL+M;
end

function [inCon,eqCon] = cons(w,x_k,t_k,f,h,R_k,T_h,N)
u_guess = reshape(w,size(R_k,1),N+1);
u = @(x,t)control(t,t_k,u_guess,T_h,N);
    function ut = control(t,t_k,u_guess,T_h,N)
        i = floor(t*N/T_h)+1;
        ut = u_guess(:,i)+(u_guess(:,i+1)-u_guess(:,i))*(t-t_k)/T_h;
    end
x_sim = dynSim(f,u,x_k,T,dt);
n_h = size(h(x_sim(:,1),u_guess(:,1)),1);
inCon = zeros(n_h*N,1);
for j = 1:N,
    inCon(:,(j-1)*n_h+(1:n_h)) = h(x_sim(:,j),u_guess(:,j));
end

eqCon = 0;
end