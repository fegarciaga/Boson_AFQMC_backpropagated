function [Phi, w, O, B, B_W, parents]=pop_cntrl_bos(Phi, w, O, N_wlk, N_sites, B, B_W, parents, N, t_bp, flag_bbp)
    %% Preparation
    % Create empty matrices that will hold the outputs
    new_Phi=zeros(2*N_sites, N_wlk); %in the end the number of walkers will still be N_wlk
    new_O=zeros(N_wlk,1);
    new_B=ones(t_bp+1, N_wlk, 2*N_sites,N+1);
    new_B_W=ones(t_bp+1, N_wlk, N+1);
    new_parents=zeros(N_wlk,N+1);
    % scaling factor to bring the current total weight back to the original level (=N_wlk)
    d=N_wlk/sum(w);
%     display(sum(w));
%     start the "comb" at a random position to avoid bias against the first walker
    sum_w=-rand;
    n_wlk=0;
    %% Apply the comb
    for i_wlk=1:N_wlk
        sum_w=sum_w+w(i_wlk)*d;
        n=ceil(sum_w);
        for j=(n_wlk+1):n
            new_Phi(:,j)=Phi(:,i_wlk);
            new_O(j)=O(i_wlk);
            % for copied orbitals, the ancestral path is actualized aswell as the corresponding father
            for i=1:N+1
                if flag_bbp(i)==1
                    new_B(:,j,:,i)=B(:,i_wlk,:,i);
                    new_B_W(:,j,i)=B_W(:,i_wlk,i);
                    new_parents(j,i)=parents(i_wlk,i);
                end
            end
        end
        n_wlk=n;
    end
    %% Return the new population, weights and overlaps:
    Phi=new_Phi;
    O=new_O;
    % All new walkers have weights to 1 and the total weight = N_wlk
    w=ones(N_wlk,1);
    B=new_B;
    B_W=new_B_W;
    parents=new_parents;
end