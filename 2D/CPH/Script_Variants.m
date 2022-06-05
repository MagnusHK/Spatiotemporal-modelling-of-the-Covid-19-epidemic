%% Script for different variants

%% Image and stencil setup
clear; close all; clc


%path = 'C:\Users\Lauri\Desktop\SimulationData\';
%path = 'C:\Users\Magnus\Desktop\DTU\Bachelorprojekt\Simulation\';
path = 'E:\Bachelor\Sim1\';


im = rgb2gray(imread("CPH_all.jpg"));


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

% Municipality images
imCPH = flip(imresize(rgb2gray(imread('Koebenhavn.jpg')), [m,n]),2) > 0;
imFred = flip(imresize(rgb2gray(imread('Frederiksberg.jpg')), [m,n]),2)>0;
imTaarn = flip(imresize(rgb2gray(imread('Taarnby.jpg')), [m,n]),2)>0;
imDrag = flip(imresize(rgb2gray(imread('Dragor.jpg')), [m,n]),2)>0;
imHvid = flip(imresize(rgb2gray(imread('Hvidovre.jpg')), [m,n]),2)>0;
imBroen = flip(imresize(rgb2gray(imread('Broendby.jpg')), [m,n]),2)>0;
imRoed = flip(imresize(rgb2gray(imread('Roedovre.jpg')), [m,n]),2)>0;
imGlo = flip(imresize(rgb2gray(imread('Glostrup.jpg')), [m,n]),2)>0;
imHer = flip(imresize(rgb2gray(imread('Herlev.jpg')), [m,n]),2)>0;
imGlad = flip(imresize(rgb2gray(imread('Gladsaxe.jpg')), [m,n]),2)>0;
imGen = flip(imresize(rgb2gray(imread('Gentofte.jpg')), [m,n]),2)>0;


% We now have a template for doing calculations, as the image is resized,
% and on a form where we can easily determine boundaries.

%% Making model parameters
Ntrue = 638117 + 103677 + 42670 + 14569 + 53451 + 35232 + ...
    41113 + 22979 + 28913 + 69200 + 74550;

gamma = 1/6;
alpha = 0;

p = [gamma; gamma; alpha];  %Ensures same gamma for I_d and I_o

%make sure that the diffusion is the image provided, with accurate values
%Doing this, we can also neglect the image as an input to functions, saving
%space
D_d = 7.5;
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
s_d = [0.6367, 0.0898, 0.0125, 0.0018, 0.0213, 0.0116, 0.0133, 0.0098, 0.0053, 0.0142, 0.0454]*1e-3;

I_d = I_d + FindICS(imCPH, s_d(1), N);  %Starting In CPH 
I_d = I_d + FindICS(imFred, s_d(2), N); %Starting in Frederiksberg
I_d = I_d + FindICS(imTaarn, s_d(3), N); %Starting in Taarnby
I_d = I_d + FindICS(imDrag, s_d(4), N);
I_d = I_d + FindICS(imHvid, s_d(5), N);
I_d = I_d + FindICS(imBroen, s_d(6), N);
I_d = I_d + FindICS(imRoed, s_d(7), N);
I_d = I_d + FindICS(imGlo, s_d(8), N);
I_d = I_d + FindICS(imHer, s_d(9), N);
I_d = I_d + FindICS(imGlad, s_d(10), N);
I_d = I_d + FindICS(imGen, s_d(11), N);
I_d = I_d .* S;

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
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\CPH');

if t_o_start == tspan(1)
    runsim(@SIR_Boundary_Variants, tspan, y0, im, beta, p, D, m, n, h, path, MakeRandom);
else
    tspan_delta = [tspan(1), t_o_start];

    %Handling new initial condition
    y0 = runsim(@SIR_Boundary_Variants, ...
                tspan_delta, y0, im, beta, p, D, m, n, h, path, MakeRandom);
    
    %%% Starting Omicron %%%
    I_o = reshape(y0(2*m*n+1:3*m*n), [m,n]);
    S = reshape(y0(1:m*n), [m,n]);

    %%%% Start of Omicron %%%%
    dose = 2/Ntrue *N; %This is a number lower than 1 (Hopefully)

    %I_o(140, 115) = dose;    %Starting Omikron in the middle at m = 264, n = 214
    I_o(259, 56) = dose;      %Starting Omicron in Dragør
    I_o = I_o(:);

    %Setting omicron into the system
    y0(1:m*n) = y0(1:m*n) - I_o;
    y0(2*m*n+1:3*m*n) = I_o;
    clear I_o S;

    %New time span
    tspan_omicron = [t_o_start, tspan(end)];

    %Running simulation for omicron
    runsim(@SIR_Boundary_Variants, ...
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

% Video = VideoWriter('DenmarkSimulation_18');
% Video.FrameRate = 30;
% open(Video)

I_d_tot = [];
I_o_tot = [];
S_tot = [];
R_tot = [];
Inew_tot = [];
t_tot = [];

%%% For municipalities %%%
S_Muni = [];
Id_Muni = [];
Io_Muni = [];
R_Muni = [];

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

    [xs, ys] = meshgrid(0:0.1105:(n-1)*0.1105, 0:0.1105:(m-1)*0.1105);

    %%% Making Data for municipalities %%%
    for i = 1:L
        Stmp = S(:,:,i);
        Idtmp = I_d(:,:,i);
        Iotmp = I_o(:,:,i);
        Rtmp = R(:,:,i);
        
        S_Muni = [S_Muni; [sum(Stmp(imCPH), 'all'), ...
                           sum(Stmp(imFred), 'all'), ...
                           sum(Stmp(imTaarn), 'all'), ...
                           sum(Stmp(imDrag), 'all'), ...
                           sum(Stmp(imHvid), 'all'), ...
                           sum(Stmp(imBroen), 'all'), ...
                           sum(Stmp(imRoed), 'all'), ...
                           sum(Stmp(imGlo), 'all'), ...
                           sum(Stmp(imHer), 'all'), ...
                           sum(Stmp(imGlad), 'all'), ...
                           sum(Stmp(imGen), 'all')]*h^2];

        Id_Muni = [Id_Muni; [sum(Idtmp(imCPH), 'all'), ...
                           sum(Idtmp(imFred), 'all'), ...
                           sum(Idtmp(imTaarn), 'all'), ...
                           sum(Idtmp(imDrag), 'all'), ...
                           sum(Idtmp(imHvid), 'all'), ...
                           sum(Idtmp(imBroen), 'all'), ...
                           sum(Idtmp(imRoed), 'all'), ...
                           sum(Idtmp(imGlo), 'all'), ...
                           sum(Idtmp(imHer), 'all'), ...
                           sum(Idtmp(imGlad), 'all'), ...
                           sum(Idtmp(imGen), 'all')]*h^2];

        Io_Muni = [Io_Muni; [sum(Iotmp(imCPH), 'all'), ...
                           sum(Iotmp(imFred), 'all'), ...
                           sum(Iotmp(imTaarn), 'all'), ...
                           sum(Iotmp(imDrag), 'all'), ...
                           sum(Iotmp(imHvid), 'all'), ...
                           sum(Iotmp(imBroen), 'all'), ...
                           sum(Iotmp(imRoed), 'all'), ...
                           sum(Iotmp(imGlo), 'all'), ...
                           sum(Iotmp(imHer), 'all'), ...
                           sum(Iotmp(imGlad), 'all'), ...
                           sum(Iotmp(imGen), 'all')]*h^2];

        R_Muni = [R_Muni; [sum(Rtmp(imCPH), 'all'), ...
                           sum(Rtmp(imFred), 'all'), ...
                           sum(Rtmp(imTaarn), 'all'), ...
                           sum(Rtmp(imDrag), 'all'), ...
                           sum(Rtmp(imHvid), 'all'), ...
                           sum(Rtmp(imBroen), 'all'), ...
                           sum(Rtmp(imRoed), 'all'), ...
                           sum(Rtmp(imGlo), 'all'), ...
                           sum(Rtmp(imHer), 'all'), ...
                           sum(Rtmp(imGlad), 'all'), ...
                           sum(Rtmp(imGen), 'all')]*h^2];
    end
    
    fprintf('Loaded file %d of %d\n', f, length(files));
    
    %%% Plotting and making a video of results %%%
%     figure(1)
%     for i = 1:L
%         Itmp = I_d(:,:,i);
%         Itmp2 = I_o(:,:,i);
%         Itmp(im>0) = Itmp(im>0) + 0.01;
%         Itmp2(im>0) = Itmp2(im>0) + 0.01;
% 
%         s = surf(xs, ys, Itmp(1:end,1:end), 'FaceAlpha',0.5);
%         hold on
%         s2 = surf(xs, ys, Itmp2(1:end,1:end), 'FaceAlpha',0.5);
%         s3 = surf(xs, ys, 0.005*ones(size(Itmp)));
%         legend('$I_\delta$', '$I_o$', Interpreter='latex')
%         xlabel('$x$ [km]', Interpreter='latex')
%         ylabel('$y$ [km]', Interpreter='latex')
%         zlabel('$I$', Interpreter='latex')
%         xlim([0, (n-1)*0.1105])
%         ylim([0, (m-1)*0.1105])
%         zlim([0, 0.5])
%         title(['SIR model with diffusion on selected municipalities for Day ' ...
%             , num2str(floor(t(i)))], Interpreter="latex")
%         view(1.925970671890411e+02,57.058355185350123)
%         s.EdgeColor = 'none';
%         s.FaceColor = [1, 0.5, 0];
%         s2.EdgeColor = 'none';
%         s2.FaceColor = [0, 1, 1];
%         s3.EdgeColor = 'none';
%         s3.FaceColor = 'b';
%         hold off
% 
%         if floor(t(i)) == 670
%             pause();
%         end
% 
%         frame = getframe(gcf);
%         writeVideo(Video, frame);
%     end

end

% close(Video)
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\CPH');

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

%% Plot of a single day

[xs, ys] = meshgrid(0:0.1105:(n-1)*0.1105, 0:0.1105:(m-1)*0.1105);

load(files(end).name)
t = data(:,1)';

day_indices = find(t == floor(t));
t = t(day_indices);
y = data(day_indices, 2:end)';

clear data;


L = length(y(1,:));
i = L;

%%% Find compartments %%%
S = reshape(y(1:m*n,:), [m,n,L]);
I_d = reshape(y(m*n+1:2*m*n,:), [m,n,L]);
I_o = reshape(y(2*m*n+1:3*m*n,:), [m,n,L]);
R = reshape(y(3*m*n+1:4*m*n,:), [m,n,L]);

clear y;

Itmp = I_d(:,:,i);
Itmp2 = I_o(:,:,i);
Itmp(im>0) = Itmp(im>0) + 0.01;
Itmp2(im>0) = Itmp2(im>0) + 0.01;


%%
close all
figure(1)
s = surf(xs, ys, Itmp(1:end,1:end), 'FaceAlpha',0.5);
hold on
s2 = surf(xs, ys, Itmp2(1:end,1:end), 'FaceAlpha',0.5);
s3 = surf(xs, ys, 0.005*ones(size(Itmp)));
legend('$I_\delta$', '$I_o$', Interpreter='latex')
xlabel('$x$ [km]', Interpreter='latex')
ylabel('$y$ [km]', Interpreter='latex')
zlabel('$I$', Interpreter='latex')
xlim([0, (n-1)*0.1105])
ylim([0, (m-1)*0.1105])
zlim([0, 0.5])
title(['SIR model with diffusion on selected municipalities for Day ' ...
    , num2str(floor(t(i)))], Interpreter="latex")
view(1.925970671890411e+02,57.058355185350123)
s.EdgeColor = 'none';
s.FaceColor = [1, 0.5, 0];
s2.EdgeColor = 'none';
s2.FaceColor = [0, 1, 1];
s3.EdgeColor = 'none';
s3.FaceColor = 'b';
print -painters
