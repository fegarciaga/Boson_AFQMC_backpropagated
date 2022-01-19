function [E_ave,E_err, Obs_bp_ave, Obs_bp_err, savedFileName]=PPMC_Bos(Lx,Ly,Lz,N_par,kx,ky,kz,U,Uab,jj,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_pc,itv_nrm,itv_Em,t_bp,t_pop,suffix)
tic; % start the  timer
bosonic_initialization; % initialize internal constants, form the trial wave function and assemble the initial population of walkers
format long;

N=floor(t_bp/N_blksteps);
flag_mea=0; %determine when a measurement should take place
flag_begin_bp=0; %determine when a back propagated measurement take place
flag_bp=0; %determine when the back propagating time ends
flag_bbp=zeros(N+1,1); %determine if back propagation is being effectuated 

E=0;
W=0;
Obs=zeros(2*N_sites,1);
Corr=zeros(N_sites,1);
W_l=0;

% Preallocate arrays:
E_blk=zeros(N_blk,1); % array to store the energy measured in every block
W_blk=zeros(N_blk,1); % array to store the total weight in every block

% Preallocate back propagation-related arrays
Phi_old=zeros(2*N_sites, N_wlk, N+1);
parents=ones(N_wlk,N+1);
B=ones(t_bp+1, N_wlk, 2*N_sites, N+1);
B_W=ones(t_bp+1, N_wlk, N+1);
% Estimate how many back propagations are possible for the total proyection time
Obs_blk=zeros(N_blk-N-1,2*N_sites);
Corr_blk=zeros(N_blk-N-1,N_sites);
W_bp=zeros(N_blk-N-1,1);

% Measurement related parameters
ind_measure=0;
ind_parents=0;
ind_end=0;
list=zeros(N+1,1);
%% Equilibration phase
for i_blk=1:N_eqblk
    for j_step=1:N_blksteps
        [Phi, Phi_old, B, B_W, parents, w, O, Obs, Corr, E, W, W_l] = bos_stepwlk(Phi, Phi_old, B, B_W, parents, N_wlk, N_sites, N, w, O, Obs, Corr, n_up, n_dn, E, W, W_l, K_old, Proj_k_half, flag_mea, flag_begin_bp, flag_bp, flag_bbp, Phi_T, N_par, t_bp, t_pop, U, Uab, fac_norm, deltau, ind_parents, ind_end, list);
        if mod(j_step,itv_pc)==0
            [Phi, w, O, B, B_W, parents]=pop_cntrl_bos(Phi, w, O, N_wlk, N_sites, B, B_W, parents, N, t_bp, flag_bbp); % population control
        end
        if mod(j_step,itv_nrm)==0
            [Phi, O, w]=bos_norm(Phi, Phi_T, N_wlk, O, w); % normalization control
        end
    end
end
%% Measurement phase    
for i_blk=1:N_blk
    for j_step=1:N_blksteps
        t=(i_blk-1)*N_blksteps+j_step;
        if mod(j_step,itv_Em)==0
            flag_mea=1;
            if i_blk<N_blk-N
                flag_begin_bp=1;
                ind_parents=mod(ind_parents+1,N+1);
                flag_bbp(ind_parents+1)=1;
            end
        else
            flag_mea=0;
            flag_begin_bp=0;
        end
        if t>(N+1)*N_blksteps
            if j_step==t_bp-N*N_blksteps
                flag_bp=1;
                ind_measure=ind_measure+1;
                ind_end=mod(ind_end+1,N+1);
            else
                flag_bp=0;
            end
        end     
        for ii=1:N+1
            if flag_bbp(ii)==1
                list(ii)=list(ii)+1;
            else
                list(ii)=0;
            end
        end
        % propagate the walkers:
        if ind_measure>0
            [Phi, Phi_old, B, B_W, parents, w, O, Obs_blk(ind_measure,:), Corr_blk(ind_measure,:), E_blk(i_blk), W_blk(i_blk), W_bp(ind_measure)] = bos_stepwlk(Phi, Phi_old, B, B_W, parents, N_wlk, N_sites, N, w, O, Obs_blk(ind_measure,:), Corr_blk(ind_measure,:), n_up, n_dn, E_blk(i_blk), W_blk(i_blk), W_bp(ind_measure), K_old, Proj_k_half, flag_mea, flag_begin_bp, flag_bp, flag_bbp, Phi_T, N_par, t_bp, t_pop, U, Uab, fac_norm, deltau, ind_parents+1, ind_end+1, list);
        else
            [Phi, Phi_old, B, B_W, parents, w, O, Obs, E_blk(i_blk), W_blk(i_blk), W] = bos_stepwlk(Phi, Phi_old, B, B_W, parents, N_wlk, N_sites, N, w, O, Obs, Corr, n_up, n_dn, E_blk(i_blk), W_blk(i_blk), W, K_old, Proj_k_half, flag_mea, flag_begin_bp, flag_bp, flag_bbp, Phi_T, N_par, t_bp, t_pop, U, Uab, fac_norm, deltau, ind_parents+1, ind_end+1, list);
        end
        if mod(j_step,itv_pc)==0
            [Phi, w, O, B, B_W, parents]=pop_cntrl_bos(Phi, w, O, N_wlk, N_sites, B, B_W, parents, N, t_bp, flag_bbp); % population control
        end
        if mod(j_step,itv_nrm)==0
            [Phi, O, w]=bos_norm(Phi, Phi_T, N_wlk, O, w); % normalization control
        end
        if mod(j_step, itv_Em)==0
            % update the exponent of the pre-factor exp(-deltau*(H-E_T))
            fac_norm=real(E_blk(i_blk)/W_blk(i_blk))*deltau+(0.5*U*real((n_up'*n_up)+(n_dn'*n_dn))+0.5*Uab*real(n_up'*n_dn+n_dn'*n_up))*deltau;
        end
        if flag_bp==1
            flag_bbp(ind_end+1)=0;    
        end
    end
    E_blk(i_blk)=E_blk(i_blk)/W_blk(i_blk);
    if ind_measure>0
        Obs_blk(ind_measure,:)=real(Obs_blk(ind_measure,:)/W_bp(ind_measure));
        Corr_blk(ind_measure,:)=real(Corr_blk(ind_measure,:)/W_bp(ind_measure));
    end
    display(strcat('E(',int2str(i_blk),')=',num2str(real(E_blk(i_blk)))))
end
%% Results
E=real(E_blk);
E_ave=mean(E)
E_err=std(E)/sqrt(N_blk)
Obs_bp_ave=mean(Obs_blk)
Obs_bp_err=std(Obs_blk)/sqrt(N_blk-N-1)
% The total computational time:
time=toc() % stops the timer

%% Save data to a *.mat file
save (savedFileName, 'E', 'E_ave', 'E_err', 'time');
save (savedFileName, '-append', 'Lx', 'Ly','Lz', 'N_par', 'kx', 'ky','kz', 'U', 'tx', 'ty','tz');
save (savedFileName, '-append', 'deltau', 'N_wlk', 'N_blksteps', 'N_eqblk', 'N_blk', 'itv_pc','itv_nrm','itv_Em');
save (savedFileName, '-append', 'K', 'Phi_T');