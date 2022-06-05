function dydt = SIR_variant(t, y, beta, gamma)
%
% beta = [beta_d, beta_o]


dydt = [-beta(1)*y(1)*y(2)-beta(2)*y(1)*y(3);
        beta(1)*y(1)*y(2)-gamma*y(2);
        beta(2)*y(1)*y(3)-gamma*y(3)];

end