function dydt = SIR_diff2(t, y, C, p, bcs)

% p(1) = beta
% p(2) = gamma
% p(3) = h


dydt = zeros(length(y), 1);
S = y(1:2:end);
I = y(2:2:end);

S_new = (C*S) - p(1)*S.*I.*bcs;
I_new = (C*I) + (p(1)*S.*I - p(2)*I).*bcs;


%Output
dydt(1:2:end) = S_new;
dydt(2:2:end) = I_new;

end