function dydt = SIR(t, y, beta, gamma)

dydt = [-beta*y(1)*y(2);beta*y(1)*y(2)-gamma*y(2)];

end