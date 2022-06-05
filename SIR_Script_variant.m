%% Simple SIR model script
clear; clc; close all

%%% Data loading %%%
table = readmatrix('Data/Municipality_cases_time_series.csv');
table = table(1:end, [3, 5, 82, 11, 43, 64, 52, 76, 39, 35, 8]);
%table = table(:,2:end);
infected = table(1:end,:); % tested positive for the first time on a given day

munis = ["Copenhagen", "Frederiksberg", "T\a rarnby", "Drag\o r", "Hvidovre", ...
          "Br\o ndby", "R\o dovre", "Glostrup", "Herlev", "Gladsaxe", "Gentofte"];

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
N = 1124471;     %Hovedstaden
%N = 5.806e6;      %Danmark

infec_o = infec(640:818, :)/N;
infec_o_tot = sum(infec_o, 2);

infec = infec(491:640, :)/N;
infec_tot = sum(infec, 2);

%% SIR model Finding an optimal beta for delta
tspan = 491:640;
I0 = infec_tot(1);
y0 = [1-I0; I0];

betalist = linspace(0, 1, 1000);
gamma = 1/6;
err = 1e10*ones(length(betalist), 1);
iter = 0;

for beta = betalist
    iter = iter + 1;
    [t, y] = ode45(@SIR, tspan, y0, [], beta, gamma);

    if y(end, 2) >= infec_tot(end)
        err(iter) = sum((y(:, 2) - infec_tot).^2)/length(t);
    end
end
ind = find(err == min(err));
beta = betalist(ind(1));

fprintf('The optimal beta for the delta variant is %.3f\n', beta);


[t, y] = ode45(@SIR, tspan, y0, [], beta, gamma);

%% Plotting
close all

figure()
plot(t, y(:,2))
hold on
plot(t, infec_tot)
grid on
legend('Simulated infected', 'Real infected',  Interpreter = 'latex')
title(['SIR models with $\beta_\delta = ', num2str(beta) ,'$ and $\gamma_\delta = 1/6$'], ...
    Interpreter = "latex" )
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter = "latex")


%% Running for the Omicron case as well
close all

tspan_d = 491:640;
tspan_o = 640:818;

I0_o = 2/N; %Omicron start
I0_d = infec_tot(1); %Delta start

y0_d = [1 - I0_d; I0_d];

beta_o = 2.29428*beta;

[td, yd] = ode45(@SIR, tspan_d, y0_d, [], beta, gamma);

y0_o = [yd(end,1)-I0_o; yd(end,2); I0_o];

[to, yo] = ode45(@SIR_variant, tspan_o, y0_o, [], [beta, beta_o], gamma);

figure()
plot([td;to], [yd(:,2); yo(:,2)])
hold on
plot(to, yo(:, 3))
plot([td;to], [infec_tot; infec_o_tot])
plot([td;to], [yd(:,2); yo(:,2)+yo(:,3)])
grid on
legend('Simulated $\delta$-variant', 'Simulated Omicron-variant', ...
        'Real infected', 'Simulated infected total',  Interpreter = 'latex')
title(['SIR model with competition of variants'], ...
    Interpreter = "latex" )
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter = "latex")


om = [zeros(size(td)); yo(:,3)];
delt = [yd(:,2); yo(:,2)];
t = [td;to];

p_o = om./(om+delt);
p_d = delt./(om+delt);

figure()
area(t,[p_o,p_d])
ylim([0,1])
xlim([490,818])
xlabel('Days', Interpreter='latex')
title('Proportion of each variant', Interpreter='latex')
legend({'Omicron', 'Delta'}, Interpreter="latex")

%%
figure()
hold on
plot(t, [yd(:,1); yo(:,1)])
plot([td;to], [yd(:,2); yo(:,2)])
plot(to, yo(:, 3))
plot([td;to], [yd(:,2); yo(:,2)+yo(:,3)])
plot(t, 1-([yd(:,2); yo(:,2)+yo(:,3)] + [yd(:,1); yo(:,1)]))
grid on
xlim([t(1), t(end)])
legend('$S$', '$I_\delta$', '$I_o$', 'Total Infections', '$R$',  Interpreter = 'latex')
title(['SIR model with competition of variants on selected municipalities'], ...
    Interpreter = "latex" )
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter = "latex")

%% Plot with diffusion (Only works when Script_variants have been run)
close all

figure()
plot(t, [yd(:,2); yo(:,2)+yo(:,3)])
hold on
%plot(t_tot, I_tot/33555)
plot(t_tot, I_tot/169165)
plot(t, [infec_tot; infec_o_tot])
grid on
xlim([491, 818])
legend(['Infected using SIR'], ...
       ['Infected using SIR with Diffusion'], ...
       'Real infected',  Interpreter = 'latex')
title(['SIR models with competition of variants'], ...
    Interpreter = "latex" )
xlabel('Days', Interpreter='latex')
ylabel('Percentage of people', Interpreter = "latex")


figure()
plot(t, abs([yd(:,2); yo(:,2)+yo(:,3)] - [infec_tot; infec_o_tot]))
hold on
%plot(t_tot, abs(I_tot'/33555 - [infec_tot; infec_o_tot]))
plot(t_tot, abs(I_tot'/169165 - [infec_tot; infec_o_tot]))
grid on
xlim([491, 818])
legend(['SIR, MSE = $8.37\cdot 10^{-4}$'], ...
       ['SIR with Diffusion, MSE = $2.33 \cdot 10^{-4}$'], ...
       Interpreter = 'latex')
xlabel('Days', Interpreter='latex')
title('Aboslute error from data at a given day', Interpreter='latex')

%% Plotting municipalities (Only works when Script_Variants have been run)
close all

for i = 1:length(Io_Muni(1,:))
    figure()
    plot(t_tot, Id_Muni(:,i)/33555)
    hold on
    grid on
    plot(t_tot, Io_Muni(:,i)/33555)
    plot(t, [infec(:,i); infec_o(:,i)])
    plot(t_tot, Id_Muni(:,i)/33555 + Io_Muni(:,i)/33555)
    ylim([0, inf])
    xlim([491, 818])
    legend('Simulated $\delta$-variant', 'Simulated omicron variant', ...
            'Real infected', 'Simulated infected total', ...
            Interpreter='latex', Location='best')
    xlabel('Days', Interpreter='latex')
    ylabel('Percentage of people', Interpreter='latex')
    title(strcat("SIR model with diffusion for ", munis(i)), Interpreter = 'latex')
end
