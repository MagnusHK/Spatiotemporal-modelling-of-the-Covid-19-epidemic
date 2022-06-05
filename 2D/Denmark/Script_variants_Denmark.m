%% Script for different variants

%% Image and stencil setup
clear; close all; clc


%path = 'C:\Users\Lauri\Desktop\SimulationData\';
%path = 'C:\Users\Magnus\Desktop\DTU\Bachelorprojekt\Simulation\';
path = 'E:\Bachelor\Sim3_Danmark\';


im = rgb2gray(imread("danmarkmotorveje.jpg"));


[r, c] = size(im);

xspan = [0, c];
yspan = [0, r];
tspan = [491, 818];
t_o_start = 640;

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

im = flip(imresize(im, [m, n]),2);

im = double(im)/255;
im(im <= 0.1) = 0;
im((im < 0.8) & (im > 0)) = 7.5;
im(im >= 0.8) = 7.5;

% image outlining cities
imby = flip(imresize(rgb2gray(imread('danmarkby.jpg')), [m,n]),2) > 230;


% We now have a template for doing calculations, as the image is resized,
% and on a form where we can easily determine boundaries.

%% Making model parameters
Ntrue = 5.806e6;

gamma = 1/6;
alpha = 0;

p = [gamma; gamma; alpha];  %Ensures same gamma for I_d and I_o

%make sure that the diffusion is the image provided, with accurate values
%Doing this, we can also neglect the image as an input to functions, saving
%space
D_d = 0.2084;
D = [D_d, 2.29428*D_d];

im = im > 0;

N = sum(im, "all");

beta_d = 0.18719/N * sum(im, 'all');

beta = [beta_d; 
        2.29428*beta_d];

MakeRandom = false; %Set to true if you wish to add random infection

clear im_cities;


%% Initial conditions

S = zeros(m, n);
S(im) = N/sum(im, "all");

I_d = zeros(m,n);
I_o = zeros(m,n);

%%% Starting values (places) at t = 491 for delta %%%
s = 3.2449e-04;

I_d = I_d + FindICS(imby, s, N);  %Starting In different cities

%clear imCPH imFred imTaarn imDrag imHvid imBroen imRoed imGlo imHer imGlad imGen

if t_o_start == tspan(1)
    I_o(150, 50) = 0.05*S(150, 50);    %Starting Omikron in bottom left ish at m = 264, n = 214
end

I_d = I_d(:);
I_o = I_o(:);

S = S(:);
R = zeros(m*n,1);

y0 = [S-I_d-I_o; I_d; I_o; R;];

clear S I_d I_o R;

%% Running simulation
cd(path);
delete *.mat
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\Denmark');

if t_o_start == tspan(1)
    runsim(@SIR_B_Variants, tspan, y0, im, beta, p, D, m, n, h, path, MakeRandom);
else
    tspan_delta = [tspan(1), t_o_start];

    %Handling new initial condition
    y0 = runsim(@SIR_B_Variants, ...
                tspan_delta, y0, im, beta, p, D, m, n, h, path, MakeRandom);
    
    %%% Starting Omicron %%%
    I_o = reshape(y0(2*m*n+1:3*m*n), [m,n]);
    S = reshape(y0(1:m*n), [m,n]);

    %%%% Start of Omicron %%%%
    dose = 2/Ntrue *N; %This is a number lower than 1 (Hopefully)

    I_o(448, 546) = dose;      %Starting Omicron in copenhagen
    I_o = I_o(:);

    %Setting omicron into the system
    y0(1:m*n) = y0(1:m*n) - I_o;
    y0(2*m*n+1:3*m*n) = I_o;
    clear I_o S;

    %New time span
    tspan_omicron = [t_o_start, tspan(end)];

    %Running simulation for omicron
    runsim(@SIR_B_Variants, ...
                tspan_omicron, y0, im, beta, p, D, m, n, h, path, MakeRandom);
end


clear y0;

%% Reformulating output and finding statistics on the problem
close all;

cd(path)
files = dir('*.mat');

T = struct2table(files); % convert the struct array to a table
T = sortrows(T, 'date'); % sort the table by 'date'
files = table2struct(T); % change it back to struct array if necessary
clear T;

%%

Video = VideoWriter('DenmarkSimulation_19');
Video.FrameRate = 30;
open(Video)

I_d_tot = [];
I_o_tot = [];
S_tot = [];
R_tot = [];
Inew_tot = [];
t_tot = [];

[xs, ys] = meshgrid((n-1)*0.6632:-0.6632:0, 0:0.6632:(m-1)*0.6632);

for f = 1:length(files)

    load(files(f).name);
    

    t = data(:,1)';
%     jump = ceil( length( data(:,1) )/( (t(end)-t(1))*16 ) );
%     
%     t = t(1:jump:end);
%     y = data(1:jump:end, 2:end)';
    day_indices = find(t == floor(t));
    t = t(day_indices);
    y = data(day_indices, 2:end)';
    
    clear data;
    
    
    L = length(y(1,:));
    
    %%% Find compartments %%%
    S = reshape(y(1:m*n,:), [m,n,L]);
    I_d = reshape(y(m*n+1:2*m*n,:), [m,n,L]);
    I_o = reshape(y(2*m*n+1:3*m*n,:), [m,n,L]);
    R = reshape(y(3*m*n+1:4*m*n,:), [m,n,L]);

    clear y;

    %%% Integration step %%%
    S_tot = [S_tot, reshape(sum(sum(S)),[1,L])*h^2];
    I_d_tot = [I_d_tot, reshape(sum(sum(I_d)),[1,L])*h^2];
    I_o_tot = [I_o_tot, reshape(sum(sum(I_o)),[1,L])*h^2];
    R_tot = [R_tot, reshape(sum(sum(R)),[1,L])*h^2];
    t_tot = [t_tot, t];

    
    fprintf('Loaded file %d of %d\n', f, length(files));
    
    %%% Plotting and making a video of results %%%
    figure(1)
    for i = 1:L
        Itmp = I_d(:,:,i);
        Itmp2 = I_o(:,:,i);
        Itmp(im>0) = Itmp(im>0) + 0.01;
        Itmp2(im>0) = Itmp2(im>0) + 0.01;

        s = surf(xs, ys, Itmp(1:end,1:end), 'FaceAlpha',0.5);
        hold on
        s2 = surf(xs, ys, Itmp2(1:end,1:end), 'FaceAlpha',0.5);
        s3 = surf(xs, ys, 0.005*ones(size(Itmp)));
        legend('$I_\delta$', '$I_o$', Interpreter='latex')
        xlabel('$x$ [km]', Interpreter='latex')
        ylabel('$y$ [km]', Interpreter='latex')
        zlabel('$I$', Interpreter='latex')
        xlim([0, (n-1)*0.6632])
        ylim([0, (m-1)*0.6632])
        zlim([0, 0.5])
        title(['SIR model with diffusion on Denmark for Day ', num2str(floor(t(i)))], Interpreter="latex")
        view(1.925970671890411e+02,57.058355185350123)
        s.EdgeColor = 'none';
        s.FaceColor = [1, 0.5, 0];
        s2.EdgeColor = 'none';
        s2.FaceColor = [0, 1, 1];
        s3.EdgeColor = 'none';
        s3.FaceColor = 'b';
        hold off
        if floor(t(i)) == 670
            pause();
        end

        frame = getframe(gcf);
        writeVideo(Video, frame);
    end

end

close(Video)
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\Denmark');

Newly = (I_d_tot + I_o_tot + R_tot)/N;
t = t_tot(2:end);
dt = t(2:end)-t(1:end-1);

Newly = (Newly(3:end)-Newly(1:end-2))./dt;

I_tot = (I_d_tot + I_o_tot);

%%
close all


figure(2)
plot(t_tot, S_tot/N)
hold on
plot(t_tot, I_d_tot/N)
plot(t_tot, I_o_tot/N)
plot(t_tot, I_tot/N)
plot(t_tot, R_tot/N)
%plot(t(1:end-1), Newly)
plot(t_tot, ones(1,length(t_tot)), 'k--')
grid on
xlim([491, 818])
legend('$S$', '$I_\delta$', '$I_o$', 'Total Infections', '$R$', Interpreter='latex')
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter='latex')
title('Competition of variants with diffusion on selected Municipalities', ...
        Interpreter='latex')


p_o = I_o_tot./(I_o_tot+I_d_tot);
p_d = I_d_tot./(I_o_tot+I_d_tot);
figure(3)
area(t_tot,[p_o',p_d'])
ylim([0,1])
xlim([491, 818])
xlabel('Days', Interpreter='latex')
title('Proportion of each variant with diffusion on selected municipalities', Interpreter='latex')
legend({'Omicron', 'Delta'}, Interpreter="latex")


