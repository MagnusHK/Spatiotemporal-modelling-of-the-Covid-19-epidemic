%% Practical convergence test of our scheme in 1D %%
clear; close all; clc;

% initialize parameters
d = 0.95; % diffusivity
a = 0;
b = pi;
m = 30;
h = (b-a)/(m+1); % mesh grid
k = 100*h.^2; % time step

% define time- and spatial domain
space = a:h:b;
time = 0:k:100;

% true solution, and plot
U=@(x,t) sin(x).*exp(-t*(d-1));
[X,T]=meshgrid(space,time);
figure(1)
surf(X,T,U(X,T));
xlabel('x')
ylabel('t')
title('True solution')


%% Defining stencil matrix A
clc;

% Initilize parameters
tspan = [0, 100];
xspan = [0, pi];

m = 30;  %Size of matrix. Number of points
a = -2; %Diagonal
b = 1;  %Superdiagonal
c = 1;  %Subdiagonal
d = 0.95; % diffusity
h = (xspan(2) - xspan(1))/(m+1);

% Stencil matrix
A = (1/h^2)*(diag(a*ones(1,m)) + diag(b*ones(1,m-1),1) + diag(c*ones(1,m-1),-1));

% Initial conditions
space = xspan(1)+h:h:xspan(2)-h;
u_0 = sin(space)';




%%% Running simulation %%%
[t, y] = ode45(@conv_test, time, u_0, [], A, d);
y_sol = zeros(size(y,1),size(y,2)+2);
y_sol(:,2:end-1) = y;

% Plot of solution
figure(2);
surf(xspan(1):h:xspan(2),t,y_sol);
xlabel('x')
ylabel('t')
title('Our solution')


%% Convergence test
clear; close all; clc;

% Initialization
err = zeros(4,1);
a = 0;
b = pi;
tspan = [0,100];
stat.k = [];
stat.h = [];
iter = 0;
d = 0.95; % diffusity
a1 = -2; %Diagonal
b1 = 1;  %Superdiagonal
c1 = 1;  %Subdiagonal
% true solution
U=@(x,t) sin(x).*exp(-t*(d-1));

for m = [10,20,40,80]
    iter = iter + 1;
    iter
    h = (b-a)/(m+1); % mesh grid
    k = 100*h.^2; % time step size
    space_bc = a:h:b;
    time = tspan(1):k:tspan(2);
    
    %Initial conditions
    space = a+h:h:b-h;
    u_0 = sin(space)';

    A = (1/h^2)*(diag(a1*ones(1,m)) + diag(b1*ones(1,m-1),1) + diag(c1*ones(1,m-1),-1));
    
    [t, y] = ode45(@conv_test, time, u_0, [], A, d);
    y_sol = zeros(size(y,1),size(y,2)+2);
    y_sol(:,2:end-1) = y;

    [X,T] = meshgrid(space_bc,time);
    
    err(iter) = max(abs(y_sol-U(X,T)),[],'all');

    stat.k = [stat.k, k];
    stat.h = [stat.h, h];
end


% Convergence plot
figure();
loglog(stat.h, 0.002.*err, '.-', 'Linewidth', 3)
hold on
plot(stat.h, stat.h, '.-')
plot(stat.h, stat.h.^2, '.-')
plot(stat.h, stat.h.^3, '.-')
grid on
title('Convergence rate of our scheme', Interpreter='latex')
xlim([stat.h(5),stat.h(1)])
xlabel('$h$', Interpreter='latex');
ylabel('Error, $||\hat{u} - u^*||_\infty$', Interpreter='latex');
legend('Our rate', '$\mathcal{O}(h)$', '$\mathcal{O}(h^2)$', '$\mathcal{O}(h^3)$', Interpreter='latex', Location='best')
    

%% Practical convergence test of our scheme 2D %%
clear; close all; clc;

% initialize parameters
d = 1; % diffusivity
ax = pi;
bx = 2*pi;
ay = pi/2;
by = 3*pi/2;
m = 10;
h = (bx-ax)/(m+1); % mesh grid
k = 100*h.^2; % time step

% Define time- and spatial domain
space_x = ax:h:bx;
space_y = ay:h:by;
time = 0:k:100;

% true solution, and plot
U=@(x,y,t) t*sin(x).*cos(y);
[X,Y]=meshgrid(space_x,space_y);
figure(1)
surf(space_x,space_y,U(space_x,space_y',1));
xlabel('x')
ylabel('y')
title('True solution for t=2')


%% Our approximation %%
clear; close all; clc;

d = 1; % diffusivity
ax = pi;
bx = 2*pi;
ay = pi/2;
by = 3*pi/2;
m = 10; % Size of matrix. Number of points
h = (bx-ax)/(m+1); % mesh grid (h_x = h_y)
k = h^2; % time step

%%% Domain %%%
space_x = ax+h:h:bx-h;
space_y = ay+h:h:by-h;

%Initial conditions
U0 = @(x,y) 0.001*sin(x).*cos(y');
U_true = @(x,y,t) t*sin(x).*cos(y');

%u_0 = zeros(m^2,1);
u_0 = U0(space_x, space_y);

% Generate stencil matrix
A = FPLaplacian(m,h);

%%% Running simulation %%%
clc; % clear command window

[t, u] = ode45(@conv_test_2D, [0.001,4], u_0(:), [], d, A);
u_tmp = reshape(u',[m,m,length(t)]);
u_sol = zeros(m+2,m+2,length(t));
u_sol(2:end-1,2:end-1,:) = u_tmp;

% Plot of solution
for i = 1:length(t)
    figure(2);
    surf(ax:h:bx,ay:h:by,u_sol(:,:,i));
    xlabel('x')
    ylabel('y')
    title('Our solution')
    pause(0.015)
end


%%% Convergence test in 2D %%%
err = [];
ms = 10:40;
for m = ms
    
    h = (bx-ax)/(m+1); % mesh grid (h_x = h_y)
    
    %%% Domain %%%
    space_x = ax+h:h:bx-h;
    space_y = ay+h:h:by-h;
    
    %%% Initial %%%
    u_0 = U0(space_x, space_y);
    A = FPLaplacian(m,h);

    %%% Time stepper %%%
    [t, u] = ode45(@conv_test_2D, [0.001,4], u_0(:), [], d, A);
    u_tmp = reshape(u',[m,m,length(t)]);
    
    err_tmp = [];
    for j = 1:length(t)
        err_tmp = [err_tmp, max(max(abs(u_tmp(:,:,j) - U_true(space_x,space_y,t(j)))))];
    end
    err = [err, max(err_tmp)];

    fprintf('Done m = %d up to %d\n', m, ms(end));

end
stat.h = (bx-ax)./(ms+1);

%%% plot of convergence rate%%%
figure();
loglog(stat.h, err, '.-', 'Linewidth', 3)
hold on
plot(stat.h, stat.h, '.-')
plot(stat.h, stat.h.^2, '.-')
plot(stat.h, stat.h.^3, '.-')
grid on
title('Convergence rate of our scheme', Interpreter='latex')
xlim([stat.h(end),stat.h(1)])
xlabel('$h$', Interpreter='latex');
ylabel('Error, $||\hat{u} - u^*||_\infty$', Interpreter='latex');
legend('Our rate', '$\mathcal{O}(h)$', '$\mathcal{O}(h^2)$', '$\mathcal{O}(h^3)$', Interpreter='latex', Location='best')





