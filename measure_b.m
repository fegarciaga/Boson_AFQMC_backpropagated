function [e] = measure_b(K, phi, Phi_T, O, N_sites, N_par, U, Uab)
    %% Calculate Green's function
    G=(Phi_T(:)').'*phi(:).';
    M=diag(G);
    K_aux=K+U/2*eye(2*N_sites,2*N_sites);
    %%  Calculate energy due to one body operators:
    one_body_energy=Phi_T'*(K_aux*phi)*N_par/O;
    %% Calculate energy due to two body operators:
    two_body_energy=U/2*(M.'*M)*N_par*(N_par-1)/O^2;
    % Calculate the green's function taking into account the off diagonal terms
    L=M(1:N_sites);
    N=M(N_sites+1:2*N_sites);
    two_body_energy=two_body_energy+Uab*(L.'*N)*N_par*(N_par-1)/O^2;
    e=real(one_body_energy+two_body_energy);
end