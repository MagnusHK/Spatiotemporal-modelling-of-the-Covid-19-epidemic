function dydt = pop(t, y, D, im, h)

[m, n] = size(im);

dydt = AmultBoundary(y,im, m,n,h,D);

end