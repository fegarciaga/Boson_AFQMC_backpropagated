function [phi_u, phi_d, Aux_force, Log_W, B_u, B_d] = V_1(phi_u, phi_d, shift, n, Uab, deltau, Aux_force, Log_W, B_u, B_d, flag_bbp, list, N)
    %% propagates the walker with the chosen auxiliary field
    exp_fact=randn+shift;
    Aux_force=Aux_force+shift*exp_fact;
    phi_u=phi_u*exp(exp_fact*sqrt(-Uab*deltau/2));
    phi_d=phi_d*exp(exp_fact*sqrt(-Uab*deltau/2));
    Log_W=Log_W+0.5*shift^2-shift*exp_fact-sqrt(-Uab*deltau/2)*n*exp_fact;
    for i=1:N+1
        if flag_bbp(i)==1
            B_u(list(i),i)=B_u(list(i),i)*exp(exp_fact*sqrt(-Uab*deltau/2));
            B_d(list(i),i)=B_d(list(i),i)*exp(exp_fact*sqrt(-Uab*deltau/2));
        end
    end
end