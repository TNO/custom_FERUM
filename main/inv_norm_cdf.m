function x = inv_norm_cdf(P)

x = erfinv(2*(P-0.5))*sqrt(2);
x = norminv(P);
