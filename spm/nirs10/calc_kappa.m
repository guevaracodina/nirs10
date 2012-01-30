function kappa = calc_kappa(B3,B3x, B3y,cov_beta_r,c)
P = cov_beta_r * (B3 * c);
Px = cov_beta_r * (B3x * c);
Py = cov_beta_r * (B3y * c);
tmp_1 = P'*P; tmp_2 = tmp_1^(-1/2); tmp_3 = tmp_2^3;
u_derx = Px*tmp_2 - (P*(P'*Px))*tmp_3;
u_dery = Py*tmp_2 - (P*(P'*Py))*tmp_3;
kappa = sqrt(abs(det([u_derx'*u_derx u_derx'*u_dery; ...
                                          u_dery'*u_derx u_dery'*u_dery])));



