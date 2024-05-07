clc; clear all; close all;

% This code solves the Cahn-Hillard equation with a vortex advection term
% implemented to simulate fluid flow. It is using a 2D Laplacian operator
% and a finite difference scheme, which is forward-Euler in time. The
% initial condition is a vector populated with random values between -1 and
% 1, with a bias towards -1 (water). 


%% PARAMETERS
x_0 = 0;         % x-coordinate of the vortex center
y_0 = 0;         % y-coordinate of the vortex center
gamma = 5;       % strength of the vortex

Nx = 40;         % Number of space steps
L = 2;           % Length of Domain [-1,1] x [-1,1]
h = L / (Nx-1);  % Space step

Nt = 3000;        % Number of timesteps
dt = 0.0000001;   % Timestep (must be very small for stability)

x = linspace(-1, 1, Nx);
y = linspace(-1, 1, Nx);
[X, Y] = meshgrid(x,y);


A = make_A(Nx);  % Construct 2D Laplacian operator with periodic b.c's

% Randomly initialize vector with values between -1 and 1, with 
% bias towards -1 

bias = 0.3;               % larger value means more biased towards -1
u_init = (rand(Nx) * 2 - 1) - bias;
u_init(u_init<(-1)) = -1; % to ensure we are between -1 and 1

u_step = u_init(:);
v_step = u_step;

k = 0.02;     % surface tension mass / time^2
tau = 1;      % density / time
delta = 0.65; % length of interface

one = ones(length(u_step),1);


% Adjust width of the vortex
r = 0.65;

% Calculate the velocity field w = < w_1, w_2 > due to the vortex (sink)
% w1 = -gamma/(2*pi) * (Y - y_0) ./ ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(3/2);
% w2 = gamma/(2*pi)  * (X - x_0) ./ ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(3/2);

% Gradients of u and v analytically
dw1_dx = -gamma * (Y - y_0) .* (X - x_0) ./ (pi * ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(5/2));
dw1_dy = gamma / (2*pi) * (1 ./ ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(3/2) - 3 * (Y - y_0).^2 ./ ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(5/2));
dw2_dx = gamma / (2*pi) * (1 ./ ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(3/2) - 3 * (X - x_0).^2 ./ ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(5/2));
dw2_dy = -gamma * (X - x_0) .* (Y - y_0) ./ (pi * ((X - x_0).^2 + (Y - y_0).^2 + r^2).^(5/2));


% Record Video
vidfile = VideoWriter('testmovievortex.mp4','MPEG-4');
set(vidfile,'FrameRate',60)
open(vidfile);

%% SOLVE SYSTEM
% Forward Euler iteration in time
for t = 1:Nt
    % advection terms in x and y direction
    adv_u = u_step .* dw1_dx(:) + v_step .* dw1_dy(:);
    adv_v = u_step .* dw2_dx(:) + v_step .* dw2_dy(:);

    u_new = u_step + (k*dt/(delta^2 * tau*h^4)) * A * ((u_step.^3-u_step) - (delta^2)*A*u_step + adv_u);
    v_new = v_step + (k*dt/(delta^2 * tau*h^4)) * A * ((v_step.^3-v_step) - (delta^2)*A*v_step + adv_v);

    u_step = u_new;
    v_step = v_new;

    % COMMENT OUT TO GET RID OF VIDEO
    if mod(t,4)==0

        % Plotting
        u_plot = reshape(u_step, Nx, Nx);
        figure(1)
        pcolor(x,y,u_plot)
        set(gcf,'color', 'w')
        pbaspect([1,1,1])
        grid off
        xlabel("$x$", Interpreter='latex', fontsize=14)
        ylabel("$y$", Interpreter='latex', fontsize=14)
        xlim([-1,1])
        ylim([-1,1])
        drawnow
        set(gca,'FontSize',20,'FontName','Times')
        
        % Record Video
        getframe(gcf);
        writeVideo(vidfile,getframe(gcf));
    end
end
close(vidfile)


%% FUNCTIONS

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


function out = plot_vortex()
% This function plots the velocity field of a vortex
% This is just for visualization purposes

% Define parameters
% vortex center at (x_0, y_0)
x_0 = 0;           
y_0 = 0;           
gamma = 10; % strength of the vortex

% Domain
x = linspace(-1, 1, 10);    
y = linspace(-1, 1, 10);    
[X, Y] = meshgrid(x, y);    

% Calculate the velocity field components
% velocity = <u,v>
u =  gamma/(2*pi) * (Y - y_0) ./ ((X - x_0).^2 + (Y - y_0).^2);
v = -gamma/(2*pi) * (X - x_0) ./ ((X - x_0).^2 + (Y - y_0).^2);

% Plot the velocity field
figure;
quiver(X, Y, u, v);
xlabel('$x$', Interpreter='latex', fontsize=18);
ylabel('$y$', Interpreter='latex', fontsize=18);
title(['$\vec{w}=\langle w_{1} ,w_{2} \rangle =\langle \frac{\gamma }{2\pi }' ...
    '\frac{y-y_{0}}{( x-x_{0})^{2} +( y-y_{0})^{2}} ,-\frac{\gamma }{2\pi }' ...
    '\frac{x-x_{0}}{( x-x_{0})^{2} +( y-y_{0})^{2}} \rangle$'], ...
    Interpreter='latex', fontsize=18);
set(gcf,'color', 'w')
pbaspect([1,1,1])
axis equal;
end