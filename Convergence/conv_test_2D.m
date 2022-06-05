function dydt = conv_test_2D(t, y, d, A)
% creates a system of coupled ODEs

f = y./t + 2*d*y;

dydt = f(:) + d*A*y;


end