function [phi, Aux_force, Log_W, B] = V_up_dn(phi, shift, n, U, deltau, Aux_force, Log_W, B, flag_bbp, list, N)
    %% propagates the walker with the chosen auxiliary field
    exp_fact=randn+shift;
    Aux_force=Aux_force+shift*exp_fact;
    phi=phi*exp(exp_fact*sqrt(-U*deltau));
    Log_W=Log_W+0.5*shift^2-shift*exp_fact-sqrt(-U*deltau)*n*exp_fact;
    for i=1:N+1
        if flag_bbp(i)==1
            B(list(i),i)=B(list(i),i)*exp(exp_fact*sqrt(-U*deltau));
        end
    end
end