function [A,B] = linearizeABSingleT(f,x_traj,u_traj)
%LINEARIZEAB linearize a given dynamic system and return state matrix and
%input matrix
%   [A,B] = LINEARIZEAB(F,XT,UT) linearizes a given dynamic F to get
%   time-dependent linearized state matrix A and input matrix B.
%   A and B are function handles in term of time t.

n = size(x_traj(1),1);
m = size(u_traj(1),1);
x = sym('x',[n,1]);
u = sym('u',[m,1]);
t = sym('t',[1,1]);


M_sym_A = jacobian(f(x,u,t),x);
M_sym_B = jacobian(f(x,u,t),u);
M_sym_A = subs(M_sym_A,x,x_traj(t));
M_sym_B = subs(M_sym_B,x,x_traj(t));

A = @(t)linearize(t,x,f,x,u,x_traj,u_traj);
B = @(t)linearize(t,u,f,x,u,x_traj,u_traj);
A = @(t)A(t);
B = @(t)B(t);

    function M = linearizeA(t,w,f,x,u,x_traj,u_traj)
        M = zeros(size(x,1),numel(w),1);
        M(:,:,1) = eval(subs(M_sym_A,[u,t],[u_traj(t),t]));
    end

    function M = linearizeB(t,w,f,x,u,x_traj,u_traj)
        M = zeros(size(x,1),numel(w),1);
        M(:,:,1) = eval(subs(M_sym_B,[u,t],[u_traj(t),t]));
    end
end