function [phi, E] =GP(phi, N_stps, K, U, Uab, N_sites, N_par, deltau)
    %% implement Gross-Pitaievskii mean field calculation for trial wavefunction
    % Use similar process to the one in AFQMC
    HGP=K;
    exp_gp=expm(-HGP*deltau);
    for i=1:N_stps
        phi=exp_gp*phi;
        O=phi(:)'*phi(:);
        for j=1:N_sites 
            [rho_up, rho_dn]=density(phi(j,:), phi(j+N_sites,:), N_par, O);
            HGP(j,j)=K(j,j)+(U/2*rho_up+Uab*rho_dn);
            HGP(j+N_sites,j+N_sites)=K(j+N_sites,j+N_sites)+(U/2*rho_dn+Uab*rho_up);
        end    
        phi(:)=phi(:)/sqrt(phi(:)'*phi(:));
        O=phi(:)'*phi(:);   
        exp_gp=expm(-HGP*deltau);
    end
    E=measure_b(K, phi, phi, O, N_sites, N_par, U, Uab);
end