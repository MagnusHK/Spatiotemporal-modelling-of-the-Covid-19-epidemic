%% Simple SIR model script
%clear; clc; close all

%%% Data loading %%%
table = readmatrix('Data/Municipality_cases_time_series.csv');
%table = table(1:end, [3, 5, 82, 11, 43, 64, 52, 76, 39, 35, 8]);
table = table(:, 2:end);    %Danmark
infected = table(1:end,:); % tested positive for the first time on a given day

% # infected on a given day
infec = zeros(size(infected));
infec(1, :) = infected(1, :);
for i = 2:height(infected)
    if i <= 6 % 1/gamma (recovery time)
        infec(i,:) = infected(i,:) + infec(i-1,:);
    else
        infec(i,:) = infec(i-1,:) + infected(i,:) - infected(i-6,:);
    end
end

%Reformulate
%N = 1124471;    %Hovedstaden
N = 5.806e6;    %Danmark
infec = infec(6:491, :)/N;
infec_tot = sum(infec, 2);

%% SIR model
tspan = 6:491;
I0 = infec_tot(1);

beta1 = 0.192;
beta2 = 0.201;

y0 = [1-I0; I0];

[t, y] = ode45(@SIR, tspan, y0, [], beta1, 1/6);
[t2, y2] = ode45(@SIR, tspan, y0, [], beta2, 1/6);

%% Plotting
figure()
plot(t, y(:,2))
hold on
plot(t2, y2(:,2))
plot(t, infec_tot)
grid on
legend(['Simulated infected with $\beta = ', num2str(beta1), '$'], ...
       ['Simulated infected with $\beta = ', num2str(beta2), '$'], ...
       'Real infected',  Interpreter = 'latex')
title(['SIR models with $\gamma = 1/6$'], ...
    Interpreter = "latex" )
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter = "latex")

%% Plot with diffusion
close all

figure()
plot(t2, y2(:,2))
hold on
%plot(t_tot, I_tot/33555)   %For municipalities
plot(t_tot, I_tot/169165)
plot(t, infec_tot)
grid on
xlim([1, 491])
legend(['Infected using SIR'], ...
       ['Infected using SIR with Diffusion'], ...
       'Real infected',  Interpreter = 'latex')
title(['SIR models with $\beta = 0.201$ and $\gamma = 1/6$'], ...
    Interpreter = "latex" )
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter = "latex")


figure()
plot(t2, abs(y2(:,2) - infec_tot))
hold on
%plot(t_tot, abs(I_tot'/33555 - infec_tot)) %For municipalities
plot(t_tot, abs(I_tot'/169165 - infec_tot))
grid on
xlim([1, 491])
legend(['SIR, MSE = $3.14\cdot 10^{-5}$'], ...
       ['SIR with Diffusion, MSE = $4.25 \cdot 10^{-7}$'], ...
       Interpreter = 'latex')
xlabel('Days', Interpreter='latex')
title('Aboslute error from data at a given day', Interpreter='latex')