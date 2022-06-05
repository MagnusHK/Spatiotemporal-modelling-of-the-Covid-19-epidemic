%% 1D-diffusion model
clear; close all; clc;

% Defining parameters
N = 300;  %Size of matrix. Number of points

% time domain and space domain
tspan = [0, 400];
xspan = [-100, 100];

%Parameters
h = (xspan(2) - xspan(1))/(N+1);
beta = 0.22;
gamma = 0.05;

%Diffusion
d = @(x) 0.22*ones(size(x));
%ddx = @(x) 0.001;
%x values for plotting
x_vals = linspace(xspan(1),xspan(2),N)';
%x values for diffusion
difference = (x_vals(2)-x_vals(1))/2;
x_vals2 = [x_vals(1)-difference
           x_vals+difference];

% Parameter vector
p = [beta; gamma; h];

% Initial conditions
infected_0 = 0.01;
ics = repmat([1;0], N, 1);
starters = [1-infected_0; infected_0];
starters = repmat(starters, 10, 1);

m = floor(length(starters)/2);
ics(N-m+1:N+m) = starters;
%ics(1:m*2) = starters;
%ics(end-m*2+1:end) = starters;

%Boundary
bcs = ones(N,1);
bcs(1) = 0;
bcs(end) = 0;

%% Defining stencil matrix
%Diagonal
a = -(d(x_vals2(1:end-1)) + d(x_vals2(2:end)));
%super and sub diagonal
b = d(x_vals2(2:end-1));

%Matrix
A = diag(a) + diag(b,-1) + diag(b,1);

%Neumann conditions
A(1, 1:3) = [-h, 0, h];
A(end, end-2:end) = [h, 0, -h];

%Scaling
A = A/(h^2);



%% Running simulation

[t, y] = ode45(@SIR_diff2, tspan, ics, [], A, p, bcs);


%% Plotting
close all;
figure(1)
subplot(3, 1, 1)
    hold on
    plot(t, y(:,1), '-');
    plot(t, y(:,2), '-');
    %plot(t, 1-(y(:,1)+y(:,2)), '-');
    legend('S', 'I', 'R');
    hold off
    
subplot(3, 1, 2)
    hold on
    plot(t, y(:,15), '-');
    plot(t, y(:,16), '-');
    %plot(t, 1-(y(:,3)+y(:,4)), '-');
    legend('S', 'I', 'R');
    hold off
    
subplot(3, 1, 3)
    hold on
    plot(t, y(:,199), '-');
    plot(t, y(:,200), '-');
    %plot(t, 1-(y(:,5)+y(:,6)), '-');
    legend('S', 'I', 'R');
    hold off
    
    
figure(2)
S_tot = sum(y(:,1:2:end),2);
I_tot = sum(y(:,2:2:end),2);
hold on
plot(t, S_tot, '-')
plot(t, I_tot, '-')
plot(t, N-(S_tot + I_tot), '-')
grid on
hold off
legend('S', 'I', 'R');


%% 3d-plots
close all;
t_incr = ceil(length(t)/150);
SI_incr = ceil(2*N/150);
if mod(SI_incr, 2)
   SI_incr = SI_incr + 1; 
end

t_ind = 1:t_incr:length(t);
S_ind = 1:SI_incr:(N*2);
I_ind = 2:SI_incr:(N*2);
x_space = linspace(xspan(1),xspan(2),length(S_ind));

i = 90;% 1:127

figure(3)
subplot(3,1,1)
sgtitle('Travelling wave of infected',Interpreter = 'latex')
plot(x_space, y(t_ind(15),I_ind))
ylim([0,1])
xlim([-100,100])
subtitle(strcat('t=',num2str(round(t(t_ind(15)))),' days'),Interpreter='latex')
subplot(3,1,2)
plot(x_space, y(t_ind(50),I_ind))
ylim([0,1])
xlim([-100,100])
ylabel('Percentage of infected',FontSize=18,Interpreter='latex')
subtitle(strcat('t=',num2str(round(t(t_ind(50)))),' days'),Interpreter='latex')
subplot(3,1,3)
plot(x_space, y(t_ind(80),I_ind))
ylim([0,1])
xlim([-100,100])
xlabel('x [km]',FontSize=18,Interpreter='latex')
subtitle(strcat('t=',num2str(round(t(t_ind(80)))),' days'),Interpreter='latex')

figure(4)
surf(x_space, t(t_ind), y(t_ind,S_ind),EdgeColor="none")
zlim([0,1])
ylabel('Days',Interpreter='latex')
xlabel('km',Interpreter='latex')
title('Susceptible',Interpreter='latex')



