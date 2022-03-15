function x = inv_student_cdf(P,nu)

x = fzero('beta_student_cdf',inv_norm_cdf(P),optimset('TolX',1e-6),P,nu);
