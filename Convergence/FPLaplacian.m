function A = FPLaplacian(m,h)
    % Creates the five point laplacian stencil matrix

    e = ones(m,1);
    S = spdiags([e -2*e e], [-1 0 1], m, m);
    I = speye(m);
    A = kron(I,S)+kron(S,I);
    A=1/(h^2)*A;

end