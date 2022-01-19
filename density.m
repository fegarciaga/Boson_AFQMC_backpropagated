function [rho_up, rho_dn]= density(phi_up, phi_dn, N_par, O)
    rho_up=N_par*diag(phi_up'*phi_up)/O;
    rho_dn=N_par*diag(phi_dn'*phi_dn)/O;
end