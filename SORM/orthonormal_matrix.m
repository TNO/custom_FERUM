function R1 = orthonormal_matrix(alpha)

nrv = length(alpha);
A = fliplr(eye(nrv)');
A(:,1) = alpha;
[Q,R] = gschmidt(A);
R1 = (fliplr(Q))';
