%% Script for Denmark (Own bounnds in general)

%% Image and stencil setup
clear; close all; clc


path = 'C:\Users\Lauri\Desktop\SimulationData\';
%path = 'C:\Users\Magnus\Desktop\DTU\Bachelorprojekt\Simulation\';
%path = 'E:\Bachelor\';


im = rgb2gray(imread("CPH_all.jpg"));

[r, c] = size(im);

xspan = [0, c];
yspan = [0, r];
tspan = [6, 491];

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



% We now have a template for doing calculations, as the image is resized,
% and on a form where we can easily determine boundaries.

%% Making model parameters
Ntrue = 638117 + 103677 + 42670 + 14569 + 53451 + 35232 + ...
    41113 + 22979 + 28913 + 69200 + 74550;

gamma = 1/6;
alpha = 0;

p = [gamma; alpha];

%make sure that the diffusion is the image provided, with accurate values
%Doing this, we can also neglect the image as an input to functions, saving
%space
D = 7.5;

im = im > 0;

N = sum(im, "all");

beta = 0.197/N * sum(im, 'all');

MakeRandom = false; %Set to true if you wish to add random infection

clear im_cities;


%% Initial conditions

S = zeros(m, n);
S(im) = N/sum(im, "all");

I = zeros(m,n);

%%% Starting conditions for alpha variant %%%
I(80, 80) = 6.225149425818896e-06*N * S(80,80);  %Starting In CPH at m = 264, n = 214
I(121, 104) = 8.893070608312709e-07*N * S(121,104); %Starting in Frederiksberg

%%% Setting IC vector %%%
I = I(:);
S = S(:);
R = zeros(m*n,1);

y0 = [S-I; I; R; I];

clear S I R;

%% Running simulation
cd(path);
delete *.mat
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\CPH');

runsim(@SIR_Boundary, tspan, y0, im, beta, p, D, m, n, h, path, MakeRandom);

clear y0;

%% Reformulating output
close all;

cd(path)
files = dir('*.mat');

T = struct2table(files); % convert the struct array to a table
T = sortrows(T, 'date'); % sort the table by 'date'
files = table2struct(T); % change it back to struct array if necessary
clear T;

%%
Video = VideoWriter('DenmarkSimulation_20');
Video.FrameRate = 30;
open(Video)

I_tot = [];
S_tot = [];
R_tot = [];
Inew_tot = [];
t_tot = [];

[xs, ys] = meshgrid(0:0.1105:(n-1)*0.1105, 0:0.1105:(m-1)*0.1105);

for f = 1:length(files)

    load(files(f).name);
    

    t = data(:,1)';
    jump = 1;%ceil( length( data(:,1) )/( (t(end)-t(1))*16 ) );
    
    
    %t = t(1:jump:end);
    %y = data(1:jump:end, 2:end)';
    
    day_indices = find(t == floor(t));
    t = t(day_indices);
    y = data(day_indices, 2:end)';

    clear data;
    
    
    L = length(y(1,:));

    S = reshape(y(1:m*n,:), [m,n,L]);
    I = reshape(y(m*n+1:2*m*n,:), [m,n,L]);
    R = reshape(y(2*m*n+1:3*m*n,:), [m,n,L]);
    Inew = reshape(y(3*m*n+1:4*m*n,:), [m,n,L]);

    clear y;

    S_tot = [S_tot, reshape(sum(sum(S)),[1,L])*h^2];
    I_tot = [I_tot, reshape(sum(sum(I)),[1,L])*h^2];
    R_tot = [R_tot, reshape(sum(sum(R)),[1,L])*h^2];
    Inew_tot = [Inew_tot, reshape(sum(sum(Inew))*h^2,[1,L])];
    t_tot = [t_tot, t];
    
    fprintf('Loaded file %d of %d\n', f, length(files));
    
    figure(1)
    for i = 1:L
        Itmp = I(:,:,i);
        Itmp(im>0) = Itmp(im>0) + 0.001;

        s = surf(xs, ys, Itmp(1:end,1:end));
        xlabel('$x$ [km]', Interpreter='latex')
        ylabel('$y$ [km]', Interpreter='latex')
        zlabel('$I$', Interpreter='latex')
        xlim([0, (n-1)*0.1105])
        ylim([0, (m-1)*0.1105])
        zlim([0, 0.02])
        title(['SIR model with diffusion on selected municipalities for Day ' ...
            , num2str(floor(t(i)))], Interpreter="latex")
        view(1.973250000000000e+02,79.258354755784055)
        s.EdgeColor = 'none';
        s.FaceColor = 'interp';

        if floor(t(i)) == 300
            pause();
        end

        frame = getframe(gcf);
        writeVideo(Video, frame);
    end

end

close(Video)
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\CPH');

Newly = (I_tot + R_tot)/N;
t = t_tot(2:end);
dt = t(2:end)-t(1:end-1);

Newly = (Newly(3:end)-Newly(1:end-2))./dt;

%%
close all


figure(2)
%plot(t_tot, S_tot/N)
hold on
plot(t_tot, I_tot/N)
%plot(t_tot, 1-(I_tot+S_tot)/N)
plot(t_tot, Inew_tot/N)
%plot(t(1:end-1), Newly)
%plot(t_tot, ones(1,length(t_tot)), 'k--')
grid on
legend('I', 'Inew', 'Newly calculated')
title('Total epidemic')


%% 
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt')

tspan = 6:491;
I0 = 1.3339605912469064e-05;

gamma = 1/6;
y0 = [1-I0; I0];

[t2, y2] = ode45(@SIR, tspan, y0, [], beta, gamma);


figure()
plot(t2, y2(:,2))
hold on
plot(t_tot, I_tot/N)
grid on
legend('Infected, $I$', 'SIR Diffusion Model', Interpreter = 'latex')
title(['SIR model with $\beta_\alpha = ' num2str(beta) '$ and $\gamma_\alpha = 1/6$'], ...
    Interpreter = "latex" )
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter = "latex")

cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\CPH');