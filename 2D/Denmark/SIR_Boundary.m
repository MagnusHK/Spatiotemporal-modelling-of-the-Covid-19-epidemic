function dydt = SIR_Boundary(t, y, im, beta, p, D, m, n, h)
% Takes vector as input. First m*n elements is S. Next m*n elements is I
% beta is the transmission rate, scalar
% p defines the parameter vector
%       p(1) = gamma
%       p(2) = alpha
% D diffusion coefficient

dydt = zeros(size(y));
mn = m*n;


%Calculating S
dydt(1:mn) = -beta.*y(1:mn).*y(mn+1:2*mn) + ...
             p(2)*y(2*mn+1:3*mn) + ...
             AmultBoundary(y(1:mn), im, m, n, h, D);

%Calculating I
dydt(mn+1:2*mn) = beta.*y(1:mn).*y(mn+1:2*mn) - p(1)*y(mn+1:2*mn) + ...
             AmultBoundary(y(mn+1:2*mn), im, m, n, h, D);

%Calculating R
dydt(2*mn+1:3*mn) = p(1)*y(mn+1:2*mn) + ...
            AmultBoundary(y(2*mn+1:3*mn), im, m, n, h, D);

%Calculating Inew
dydt(3*mn+1:4*mn) = beta.*(dydt(1:mn).*y(mn+1:2*mn) + y(1:mn).*dydt(mn+1:2*mn));




end