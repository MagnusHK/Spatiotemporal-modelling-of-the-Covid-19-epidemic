%% Script for Denmark (Own bounnds in general) with population size bigger in cities

%% Image and stencil setup
clear; close all; clc


path = 'C:\Users\Lauri\Desktop\SimulationDataVeje\Sim3\';
%path = 'C:\Users\Magnus\Desktop\DTU\Bachelorprojekt\Simulation\';
%path = 'E:\Bachelor';


im = rgb2gray(imread("danmarkmotorveje.jpg"));
im_cities = rgb2gray(imread("danmarkby.jpg"));

% if length(size(im)) == 3
%     im = im(:,:,1) + im(:,:,2) + im(:,:,3) > 0;
% end

[r, c] = size(im);

xspan = [0, c];
yspan = [0, r];
tspan = [0, 500];

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
im_cities = imresize(im_cities, [m, n]);
%im = MakeInterior(m, n, im);

im = double(im)/255;
im(im <= 0.1) = 0;
im((im < 0.8) & (im > 0)) = 0.4;
im(im >= 0.8) = 1;

im_cities = double(im_cities)/255;
im_cities(im_cities <= 0.1) = 0;
im_cities((im_cities < 0.9) & (im_cities > 0)) = 0.3;
im_cities(im_cities >= 0.9) = 0.5;



% We now have a template for doing calculations, as the image is resized,
% and on a form where we can easily determine boundaries.

%% Making model parameters
N = 6e6;

gamma = 0.125;
alpha = 0.005;

p = [gamma; alpha];

%make sure that the diffusion is the image provided, with accurate values
%Doing this, we can also neglect the image as an input to functions, saving
%space
D = im;
beta = 0.4; %im_cities(:);

im = im > 0;



%% Initial conditions

%%%%%%%%%%%%%%%%%%%%%%
%%%   Suceptible   %%%
%%%%%%%%%%%%%%%%%%%%%%
% Say that 50% of our population is habitated in cities.
% For a start, they are done so uniformly, such that each point in a city
% contains the same amount of people

S = zeros(m, n);
%S(im) = N/sum(im, "all");

S(im_cities == 0.5) = (N/2)/sum(im_cities == 0.5, "all");
S(im_cities == 0.3) = (N/2)/sum(im_cities == 0.3, "all");
scale = max(S(:));

S = S(:)/scale;
pop = S;

clear im_cities;


%%%%%%%%%%%%%%%%%%%%%%
%%%    Infected    %%%
%%%%%%%%%%%%%%%%%%%%%%

I = zeros(m,n);

I(450, 69) = 0.01;    %Starting in Copenhagen 697
I(343, 347) = 0.01;    %Starting in Aarhus 697

% I(647, 97) = 0.01;    %Starting close to Copenhagen 1000
% I(490, 495) = 0.01;    %Starting close to Aarhus 1000

I = I(:);

y0 = [S-I.*S; I.*S];

clear S I;

%% Running simulation
cd(path);
delete *.mat
cd('C:\Users\Lauri\MATLAB Drive\MagLau\ML BCS projekt\2D\Denmark');

runsim(@SIR_Boundary, tspan, y0, im, beta, p, D, m, n, h, path, pop);

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
Video = VideoWriter('DenmarkSimulation_3');
Video.FrameRate = 60;
open(Video)

t_tot = [];
I_tot = [];

for f = 1:length(files)

    load(files(f).name);
    
    t = data(:,1)';
    y = data(:, 2:end)';

    t_tot = [t_tot, t];
    
    clear data;
    
    
    L = length(y(1,:));

    S = reshape(y(1:m*n,:), [m,n,L]);
    I = reshape(y(m*n+1:2*m*n,:), [m,n,L]);
    
    waitfactor = 8/L;
    
    clear y;
    
    for i = 1:L
        
        Itmp = I(:,:,i)*scale;
        I_tot = [I_tot, sum(Itmp, "all")];

        Itmp(im>0) = Itmp(im>0) + 0.3;

        figure(1)
        s = surf(Itmp(1:end,1:end));
        xlabel('Space x')
        ylabel('space y')
        %zlim([0, 1])
        view(1.973250000000000e+02,79.258354755784055)
        s.EdgeColor = 'none';
        s.FaceColor = 'interp';
        pause(waitfactor);

        frame = getframe(gcf);
        writeVideo(Video, frame);

    end

end

close(Video)

figure(2)
plot(t_tot, I_tot)
title('Total infected over time')
xlabel('Time t [days]')
ylabel('#Infected people')

