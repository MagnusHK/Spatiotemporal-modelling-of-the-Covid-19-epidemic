%% Population change script
clear; clc; close all

im = rgb2gray(imread("CPH_all.jpg"));

[r, c] = size(im);

xspan = [0, c];
yspan = [0, r];
tspan = [0, 100];

%This will put m points in the longest direction, and n points in the
%shortest direction
mdesired = 100;
m = max([r, c, mdesired]);
h = max(r,c)/(m+1);
n = length(h:h:min(r,c)-h);

% Resizing image
if max(r,c) == c
    m_temp = m;
    m = n;
    n = m_temp;
end

im = imresize(im, [m, n]);
%im = MakeInterior(m, n, im);

im = double(im)/255;
im(im <= 0.1) = 0;
im((im < 0.8) & (im > 0)) = 15;
im(im >= 0.8) = 15;


%% Parameters

rho0 = @(x, y) (sin(0.05*x).*cos(0.05*y') + 1)/2;

D = 7.5;

im = im > 0;

x = 0+h:h:n;
y = 0+h:h:m-h;

u0 = rho0(x,y);
u0(im == 0) = 0;

%% Running simulation

[t, u] = ode45(@pop, tspan, u0(:), [], D, im, h);

%%
u = u';

%% Plotting

utmp = reshape(u(:,end), [m,n]);
figure(1)
s = surf(flip(x),y,utmp);
xlabel('$x$ [pixels]', Interpreter='latex')
ylabel('$y$ [pixels]', Interpreter='latex')
title('Population distribution at day 100', Interpreter='latex')
zlim([0, 2])
xlim([min(x),max(x)])
ylim([min(y),max(y)])
s.EdgeColor = 'none';
s.FaceColor = 'Interp';
view(1.973250000000000e+02,79.258354755784055)

figure(2)
s = surf(flip(x),y,u0);
xlabel('$x$ [pixels]', Interpreter='latex')
ylabel('$y$ [pixels]', Interpreter='latex')
title('Population distribution at day 0', Interpreter='latex')
zlim([0, 2])
xlim([min(x),max(x)])
ylim([min(y),max(y)])
s.EdgeColor = 'none';
s.FaceColor = 'Interp';
view(1.973250000000000e+02,79.258354755784055)