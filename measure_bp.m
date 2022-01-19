function [O, Corr]=measure_bp(Phi_T, phi_old, B, Proj_k_half, t_bp, t_pop, N_par, N_sites)
  %% Implement backwards propagation with the stored values
  % calculate backpropagated bra
  Phi_left=Phi_T;
    for ii=1:t_bp+1
        Phi_left=Proj_k_half*Phi_left;
        for jj=1:N_sites
            Phi_left(jj)=Phi_left(jj)*B(t_bp-ii+2,jj);
            Phi_left(jj+N_sites)=Phi_left(jj+N_sites)*B(t_bp-ii+2,jj+N_sites);
        end
        Phi_left=Proj_k_half*Phi_left;
        if mod(ii,t_pop)==0
            Phi_left=Phi_left/sqrt(Phi_left'*Phi_left);
        end
    end
    %% Once the bra is back propagated, the observable is calculated the usual way
    % calculate inverse matrices in order to calculate Green's function
    G=(Phi_left(:)').'*phi_old(:).';
    O_prov=Phi_left'*phi_old;
    O=diag(G)*N_par/O_prov;
    % Pcik reference site for correlation calculations 
    i_ref=N_sites/2;
    Corr=zeros(n_sites,1);
    for i=1:N_sites
        if i==i_ref
            Corr(i)=(G(i_ref,i_ref)*G(i,i)-G(i_ref+N_sites,i_ref+N_sites)*G(i,i)-G(i_ref,i_ref)*G(i+N_sites,i+N_sites)+G(i_ref+N_Sites,i_ref+N_sites)*G(i+N_sites,i+N_sites))*N_par*(N_par-1)/O_prov^2+(G(i,i)+G(i+N_sites,i+N_sites))*N_par/O_prov;
        else
            Corr(i)=(G(i_ref,i_ref)*G(i,i)-G(i_ref+N_sites,i_ref+N_sites)*G(i,i)-G(i_ref,i_ref)*G(i+N_sites,i+N_sites)+G(i_ref+N_Sites,i_ref+N_sites)*G(i+N_sites,i+N_sites))*N_par*(N_par-1)/O_prov^2;
        end
    end
end