function dydt = conv_test(t, y, A, d)
% creates a system of coupled ODEs

dydt = y + d.*A*y;

end