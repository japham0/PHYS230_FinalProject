clc; clear all; close all;

% This code shows that our laplacian operator is accurate and
% shows that it solves the 2D heat equation with O(h^2)convergence
% in space, by solving the 2D heat equation u_t = u_xx + u_yy
% with periodic boundary conditions on the domain [-1,1] x [-1,1]
% It produces a convergence plot where the loglog slope of the error
% is parallel to the loglog plot of the space step. Since the lines are 
% parallel and both have a slope of 2, our operator has order 2 
% convergence in space. One can also see that halving the space step h
% decreases the error by a factor of 4 as expected (err vector)

Nx_list = [5, 10, 20];
h_list = [];
iter = 1;
err = [];  % Store error

while iter < 4
Nx = Nx_list(iter);
L = 2;        % Domain [-1,1] x [-1,1]
h = L / (Nx-1);
h_list(iter) = h;
Nt = 500;     % Number of timesteps
dt = 0.0000001;  
mu = 1;

x = linspace(-1, 1, Nx);
y = linspace(-1, 1, Nx);
[X, Y] = meshgrid(x,y);

u_true = zeros(Nx);
u_true_time = [];

for t = 1:Nt
    u_true = exp(-2*t*dt*pi^2)*(sin(pi*X).*sin(pi*Y));
    u_true_time = [u_true_time, u_true(:)];
end
% figure(1)
% pcolor(x,y,u_true)

A = make_A(Nx);
u_init = sin(pi*X).*sin(pi*Y);
u_step = u_init(:);

for t = 1:Nt
    u_new = u_step + (dt*mu/h^2) * A * u_step;
    u_step = u_new;
end

u_end = reshape(u_step, Nx, Nx);
u_true_end = u_true_time(:, end);

% SURFACE PLOT
% figure(2)
% surf(x,y, reshape(u_true_end, Nx, Nx))
% figure(3)
% surf(x,y, u_end - reshape(u_true_end, Nx, Nx))
% xlabel("x")
% ylabel("y")

err(iter) = norm(u_step - u_true_end, "inf")
iter = iter + 1;
end

% PLOTTING
loglog(h_list, h_list.^2)
hold on
shifting_constant = 900;
loglog(h_list, shifting_constant*err, Marker="x", MarkerEdgeColor='k', MarkerSize=14)
grid on
xlabel("h", Interpreter="latex", Fontsize=18)
ylabel("$\Vert u_{approx} - u_{true} \Vert_{\infty}$", Interpreter="latex", Fontsize=18)
legend("h", "error")

function out = make_A(Nx)
d = ones(Nx^2, 1);
A = spdiags([d d -4*d d d],[-Nx -1 0 1 Nx],...
    (Nx^2),(Nx^2));
A = full(A);
for i = 1:Nx
    for j = 1:Nx
        n = i + (j-1) * Nx;
        % RED boundary
        if i==1 && j==1
            A(n,n+(Nx-2)) = 1;
            A(n,n+(Nx^2-2*Nx))=1;
        end
        if i==Nx && j==1
            A(n,n-(Nx-2)) = 1;
            A(n,n+1) = 0;
            A(n,n+(Nx^2-2*Nx)) = 1;
        end
        if i== 1 && j== Nx
            A(n,n+(Nx-2)) = 1;
            A(n,n-1) = 0;
            A(n,n-(Nx^2-2*Nx)) = 1;
        end
        if i== Nx && j== Nx
            A(n,n-(Nx-2)) = 1;
            A(n, n-(Nx^2-2*Nx)) = 1;
        end
    end
end

for j = 2:Nx-1
    n = 1 + (j-1)* Nx;
    A(n,n+(Nx-2)) = 1;
    A(n, n-1) = 0;
end
for j = 2:Nx-1
    n = Nx*j;
    A(n,n-(Nx-2)) = 1;
    A(n, n+1) = 0;
end
for i = 2:Nx-1
    A(i,i+(Nx^2-2*Nx))=1;

    n = i + (Nx-1)*Nx;
    A(n,n-(Nx^2-2*Nx)) = 1;
end
out = A;
end

