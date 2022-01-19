function [phi, phi_old, B, B_W, parents, w, O, Obs_bp, Corr_bp, E, W, W_bp] = bos_stepwlk(phi, phi_old, B, B_W, parents, N_wlk, N_sites, N, w, O, Obs_bp, Corr_bp, n_up, n_dn, E, W, W_bp, K, Proj_k_half, flag_mea, flag_begin_bp, flag_bp, flag_bbp, Phi_T, N_par, t_bp, t_pop, U, Uab, fac_norm, deltau, ind_parents, ind_end, list)
    %% Propagate each walker:
    e=zeros(N_wlk,1); % Array containing the energy of each walker
    obs_bp=zeros(N_wlk,2*N_sites); % For the first example let's consider middle site Delta 
    corr_bp=zeros(N_wlk,N_sites);
    % stores population at a given time for future back propagation measurement
    if flag_begin_bp==1
        phi_old(:,:,ind_parents)=phi;
    end
    for i_wlk=1:N_wlk
        Phi=phi(:,i_wlk);
        if flag_begin_bp==1
            parents(i_wlk,ind_parents)=i_wlk;
        end
        if w(i_wlk)>0
            phase=N_par*angle(O(i_wlk));
            if cos(phase)<=0
                w(i_wlk)=0;
            end
            if w(i_wlk)>0
                n_up_t=N_par*diag((Phi_T(1:N_sites)').'*Phi(1:N_sites).')/O(i_wlk);
                n_dn_t=N_par*diag((Phi_T(1+N_sites:2*N_sites)').'*Phi(1+N_sites:2*N_sites).')/O(i_wlk);
                w(i_wlk)=w(i_wlk)*cos(phase);
                Aux_force=0;
                % propagate by the kinetic term exp(-1/2*deltau*K)
                [Phi]=bos_halfK(Phi, Proj_k_half);
                % propagate each lattice site of a walker by the potential term:
                Log_W=-phase*sqrt(-1);
                for j_site=1:N_sites
                    [Phi(j_site), Aux_force, Log_W, B(:,i_wlk,j_site,:)]=V_up_dn(Phi(j_site), sqrt(-1*U*deltau)*(n_up_t(j_site)-n_up(j_site)), n_up(j_site), U, deltau, Aux_force, Log_W, B(:,i_wlk,j_site,:), flag_bbp, list, N);
                end
    %                 
                for j_site=1:N_sites
                    [Phi(j_site+N_sites), Aux_force, Log_W, B(:,i_wlk,j_site+N_sites,:)]=V_up_dn(Phi(j_site+N_sites), sqrt(-1*U*deltau)*(n_dn_t(j_site)-n_dn(j_site)), n_dn(j_site), U, deltau, Aux_force, Log_W, B(:,i_wlk,j_site+N_sites,:), flag_bbp, list, N);
                end
    %                 
                for j_site=1:N_sites
                    [Phi(j_site), Phi(j_site+N_sites), Aux_force, Log_W, B(:,i_wlk,j_site,:), B(:,i_wlk,j_site+N_sites,:)]=V_1(Phi(j_site), Phi(j_site+N_sites), sqrt(-1*Uab*deltau/2)*((n_up_t(j_site)+n_dn_t(j_site))-(n_up(j_site)+n_dn(j_site))), n_up(j_site)+n_dn(j_site), Uab, deltau, Aux_force, Log_W, B(:,i_wlk,j_site,:), B(:,i_wlk,j_site+N_sites,:), flag_bbp, list, N);
                end
    %                 
                for j_site=1:N_sites
                    [Phi(j_site), Phi(j_site+N_sites), Log_W, B(:,i_wlk,j_site,:), B(:,i_wlk,j_site+N_sites,:)]=V_2(Phi(j_site), Phi(j_site+N_sites), sqrt(Uab*deltau/2)*((n_up_t(j_site)-n_dn_t(j_site))-(n_up(j_site)-n_dn(j_site))), n_up(j_site)-n_dn(j_site), Uab, deltau, Log_W, B(:,i_wlk,j_site,:), B(:,i_wlk,j_site+N_sites,:), flag_bbp, list, N);
                end
                % propagate by the kinetic term exp(-1/2*deltau*K)
                [Phi]=bos_halfK(Phi, Proj_k_half); 
                % applies the phaseless condition 
                O_new=Phi_T'*Phi;
                w(i_wlk)=w(i_wlk)*exp(fac_norm);
                x_path=max(0,cos(imag(Aux_force)));
                for i=1:N+1
                    if flag_bbp(i)==1
                        B_W(list(i),i_wlk,i)=x_path;
                    end
                end
                Log_W=Log_W+N_par*log(O_new/O(i_wlk));
                w(i_wlk)=w(i_wlk)*abs(exp(Log_W))*x_path;
                O(i_wlk)=O_new;
                if w(i_wlk)>0
                    % measure the energy if needed:
                    if flag_mea==1
                        [e(i_wlk)]=measure_b(K, Phi(:), Phi_T, O_new, N_sites, N_par, U, Uab);
                    end
                    if flag_bp==1
                        % Selects right father
                        i_father=parents(i_wlk, ind_end);
                        [obs_bp(i_wlk,:), corr_bp(i_wlk,:)]=measure_bp(Phi_T, phi_old(:,i_father,ind_end), B(:,i_wlk,:, ind_end), Proj_k_half, t_bp, t_pop, N_par, N_sites);
                    end
                end   
            end
        end
        phi(:,i_wlk)=Phi;
    end
    %% Compute the ensemble's total energy and weight if measurement took place
    if flag_mea==1
        for i_wlk=1:N_wlk
            if w(i_wlk)>0
                E=E+e(i_wlk)*w(i_wlk);
                W=W+w(i_wlk);
            end
        end
    end

    %% Compute the ensemble's given observable if measurement took place
    if flag_bp==1
        for i_wlk=1:N_wlk
            if w(i_wlk)>0
                w_aux=w(i_wlk);
                for ii=1:t_bp+1
                    w_aux=w_aux/B_W(ii,i_wlk,ind_end);
                end
                Obs_bp=Obs_bp+obs_bp(i_wlk,:)*w_aux;
                Corr_bp=Corr_bp+obs_bp(i_wlk,:)*w_aux;
                W_bp=W_bp+w_aux;
            end
        end
        B(:,:,:,ind_end)=1;
        phi_old(:,:,ind_end)=0;
        B_W(:,:,ind_end)=1;
    end
    
end