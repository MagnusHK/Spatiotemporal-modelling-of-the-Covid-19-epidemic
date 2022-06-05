%% Data script: script that utilizes the provided data %%
clear; clc; close all


%%% Data loading %%%
table = readmatrix('Data/Municipality_cases_time_series.csv');
table = table(1:end, [3, 5, 82, 11, 43, 64, 52, 76, 39, 35, 8]);
%table = table(:, 2:end);    %Danmark
infected = table(1:end,:); % tested positive for the first time on a given day

% # infected on a given day (gamma = 1/6, i.e. infected for 6 days)
infec = zeros(size(infected));
infec(1, :) = infected(1, :);
for i = 2:height(infected)
    if i <= 6 % 1/gamma (recovery time)
        infec(i,:) = infected(i,:) + infec(i-1,:);
    else
        infec(i,:) = infec(i-1,:) + infected(i,:) - infected(i-6,:);
    end
end


% Reformulate to density by using total population
N = 1124471;    %Hovedstaden
%N = 5.806e6;    %Danmark
infec = infec/N;
infec_tot = sum(infec, 2);

% plot of number of infected for a given day
figure()
plot(-1, -1)
hold on
plot(-1, -2)
plot(infec_tot)
grid on
xlim([1, length(infec_tot)])
ylim([0, inf])
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter='latex')
title('Number of infected at a certain day for selected municipalities', Interpreter='latex')



%% Data
clc; clear; close all;

% Read data matrix
table = readmatrix('Data/Municipality_cases_time_series.csv');

% fra 26-02-2020 til 24-03-2022
infected = sum(table(:,2:end), 2); % tested positive for the first time on a given day
tested = table(1:end-2,6); % tested on a given day




% # infected on a given day
infec = zeros(length(infected),1);
infec(1) = infected(1);
for i = 2:length(infected)
    if i <= 6 % 1/gamma (recovery time)
        infec(i) = infected(i) + infec(i-1);
    else
        infec(i) = infec(i-1) + infected(i) - infected(i-6);
    end
end

% infec in [0,1]
N = 5806000; % total population in Denmark
infec = infec./N;

% choosing an estimation period, here 1/7/2021 - 27/11/2021
infec = infec(491:640); % 1/7/2021 - 27/11/2021



%%% Parameter estimation %%%
% init vals
tspan = linspace(0,length(infec),length(infec));
infected_0 = infec(1);
gamma = 1/6;

% list of beta values
betalist = linspace(0,1,100);
SSE = zeros(1,length(betalist));
for i = 1:length(betalist)
    % model computations
    [t,yhat] = ode45(@SIR, tspan, [1-infected_0, infected_0],[],betalist(i), gamma);
    r = infec - yhat(:,2); % residuals
    SSE(i) = norm(r,2)^2; % sse
end


% simulation for given estimation of beta and plots
[~,ind] = min(SSE);
[t,yhat] = ode45(@SIR, tspan, [1-infected_0, infected_0],[],betalist(ind),gamma);

figure();
plot(t,yhat(:,1),t,yhat(:,2),t,1-yhat(:,1)-yhat(:,2),t,infec)
xlim([0,length(infec)])
ylim([0,1])
xlabel('t')
ylabel('Infected')
legend('S','I','R','real data')

figure()
plot(t,yhat(:,2),t,infec)
xlim([0,length(infec)])
xlabel('t')
ylabel('Infected')
legend('I','real data')
%% variant parameter estimation
clc; clear; close all;

% read data
table = readmatrix('SSI_variant_data.xlsx');

%[date, new_pos, prev_pos, tested]
% fra 27-10-2021 til 19-12-2021
omikron = table(1:end,6); % omicron infected
varianter = table(1:end,4); % other variants
varianter(end) = 0;

% # omicron infected on a given day
infec_omikron = zeros(length(omikron),1);
infec_omikron(1) = omikron(1);
for i = 2:length(omikron)
    if i <= 6 % 1/gamma (recovery time)
        infec_omikron(i) = omikron(i) + infec_omikron(i-1);
    else
        infec_omikron(i) = infec_omikron(i-1) + omikron(i) - omikron(i-6);
    end
end

% # other variant infected on a given day
infec = zeros(length(varianter),1);
infec(1) = varianter(1);
for i = 2:length(varianter)
    if i <= 6 % 1/gamma (recovery time)
        infec(i) = varianter(i) + infec(i-1);
    else
        infec(i) = infec(i-1) + varianter(i) - varianter(i-6);
    end
end

% plot
figure();
hold on
plot(infec_omikron)
plot(infec)
xlabel('t')
ylabel('Infected')
legend('Omicron infected','Other variants')
xlim([0,length(infec)])

% rescaling
% infec in [0,1]
N = 5806000/100; % nedskalering??
infec_omikron = infec_omikron./N;
infec = infec./N;


% Parameter estimation
% init vals
tspan = linspace(0,length(infec_omikron),length(infec_omikron));
infected_0 = 0.000001;

betalist = linspace(0,1,100);
SSE = zeros(1,length(betalist));
for i = 1:length(betalist)
    % model computations
    [t,yhat] = ode45(@SIRS, tspan, [1-infected_0, infected_0],[],betalist(i));
    r = infec_omikron - yhat(:,2); % residuals
    SSE(i) = norm(r,2)^2; % sse
end


%
[~,ind] = min(SSE);
[t,yhat] = ode45(@SIRS, tspan, [1-infected_0, infected_0],[],betalist(ind));

figure();
plot(t,yhat(:,1),t,yhat(:,2),t,1-yhat(:,1)-yhat(:,2),t,infec_omikron)
xlim([0,length(infec_omikron)])
ylim([0,1])
xlabel('t')
ylabel('Infected')
legend('S','I','R','real data')


%% Estimering af parametre for to varianter (delta og omikron BA.1) i én population
clc
clear
close all

table = readmatrix('SSI_variant_data.xlsx');

%[date, new_pos, prev_pos, tested]
% fra 27-10-2021 til 19-12-2021
omikron = table(1:end,6); % omicron infected
varianter = table(1:end,4); % other variants
varianter(end) = 0;


% # omicron infected on a given day
infec_omikron = zeros(length(omikron),1);
infec_omikron(1) = omikron(1);
for i = 2:length(omikron)
    if i <= 8 % 1/gamma (recovery time)
        infec_omikron(i) = omikron(i) + infec_omikron(i-1);
    else
        infec_omikron(i) = infec_omikron(i-1) + omikron(i) - omikron(i-8);
    end
end

% # other variant infected on a given day
infec = zeros(length(varianter),1);
infec(1) = varianter(1);
for i = 2:length(varianter)
    if i <= 8 % 1/gamma (recovery time)
        infec(i) = varianter(i) + infec(i-1);
    else
        infec(i) = infec(i-1) + varianter(i) - varianter(i-8);
    end
end

% infec in [0,1]
N = 5806000/100; % nedskalering??
infec_omikron = (infec_omikron./N)/2;
infec = (infec./N)/2;

% Parameter estimation
% init vals
tspan = linspace(0,length(infec_omikron),length(infec_omikron));
infected_o = infec_omikron(1);
infected_d = infec(1);
y0 = [1-infected_o-infected_d, infected_o, infected_d];

betalist = linspace(0,1,100);
SSE = zeros(length(betalist),length(betalist));

for i = 1:length(betalist)
    for j = 1:length(betalist)
        % model computations
        [~,yhat] = ode45(@SIRVariant_varbeta, tspan, y0, [], betalist(i), betalist(j));
        r_o = infec_omikron - yhat(:,2); % residuals
        r_d = infec - yhat(:,3); % residuals
        SSE(i,j) = norm([r_o,r_d],2)^2; % sse 
    end
end

%%
[i,j] = find(SSE==min(min(SSE)));
%[i1,j1] = find(SSE_d == min(min(SSE_d)));
[t,yhat] = ode45(@SIRVariant_varbeta, tspan, [1-infected_o-infected_d, infected_o, infected_d],[],betalist(i),betalist(j));

figure();
plot(t,yhat(:,1),t,yhat(:,2),t,yhat(:,3),t,infec_omikron,t,infec)
xlim([0,length(infec_omikron)])
ylim([0,1])
xlabel('t')
ylabel('Infected')
legend('S','I_o','I_\delta','I_o data','I_\delta data')





