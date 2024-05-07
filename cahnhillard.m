clc; clear all; close all;

% This code solves the Cahn-Hillard equation using a 2D Laplacian operator
% and a finite difference scheme, which is forward-Euler in time. The
% initial condition is a vector populated with random values between -1 and
% 1, with a bias towards -1 (water). 

Nx = 20;     % Grid dimensions
L = 2;       % Domain [-1,1] x [-1,1]
h = L / (Nx-1); % Space step

Nt = 100000;     % Number of timesteps
dt = 0.000001;   % Timestep (must be very small for stability)


x = linspace(-1, 1, Nx);
y = linspace(-1, 1, Nx);
[X, Y] = meshgrid(x,y);


A = make_A(Nx);  % Construct 2D Laplacian operator with periodic b.c's

% Randomly initialize vector with values between -1 and 1, with 
% bias towards -1 

bias = 0.25; % larger value of bias means more biased towards -1
u_init = (rand(Nx) * 2 - 1) - bias; 

u_step = u_init(:);

k = 0.01;  % surface tension (what units?)
tau = 0.2;   % scaling factor  (what units?)
delta = 0.1; % length of interface  (what units?)

one = ones(length(u_step),1);

% Forward Euler interation in time
for t = 1:Nt
    %u_new = u_step + (k*dt/tau)*(1/h^4)*A*(A*u_step - (1/delta^2)*u_step.*(1-u_step.^2));
    u_new = u_step + (k*dt/(delta^2 * tau*h^4)) * A * ((u_step.^3-u_step) - (delta^2)*A*u_step);
    u_step = u_new;
end

% Plotting
u_end = reshape(u_step, Nx, Nx);
figure(1)
pcolor(x,y,u_end)
xlabel("x")
ylabel("y")


% Construct the 2D Laplacian with periodic boundary conditions
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

