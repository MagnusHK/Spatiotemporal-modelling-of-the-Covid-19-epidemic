function dydt = SIR_B_Variants(t, y, im, beta, p, D, m, n, h)
%Takes vector as input. First m*n elements is S. Next m*n elements is I_d
% beta defines the vector
%       beta(1) = beta_delta
%       beta(2) = beta_o
% p defines the parameter vector
%       p(1) = gamma_delta
%       p(2) = gamma_o
%       p(2) = alpha
% D is the diffusion vector
%       D(1) = D_delta
%       D(2) = D_o

%%% Define parameters %%%
mn = m*n;
% Compartments
S = y(1:mn);
I_d = y(mn+1:2*mn);
I_o = y(2*mn+1:3*mn);
R = y(3*mn+1:4*mn);

clear y;

%Calculating S
dS = -(beta(1)*I_d + beta(2)*I_o).*S + p(3)*R + ...
            AmultBoundary(S, im, m, n, h, D(1));

%Calculating I_delta
dI_d = beta(1)*I_d.*S - p(1)*I_d + ...
            AmultBoundary(I_d, im, m, n, h, D(1));


%Calculating I_o
dI_o = beta(2)*I_o.*S - p(2)*I_o + ...
            AmultBoundary(I_o, im, m, n, h, D(2));

%Calculating R
dR = p(1)*I_d + p(2)*I_o - p(3)*R + ...
            AmultBoundary(R, im, m, n, h, D(1));


%%% Output %%%
dydt = [dS; dI_d; dI_o; dR];

end