function [phi_u, phi_d, Log_W, B_u, B_d] = V_2(phi_u, phi_d, shift, n, Uab, deltau, Log_W, B_u, B_d, flag_bbp, list, N)
    %% propagates the walker with the chosen auxiliary field
    exp_fact=randn+shift;
    phi_u=phi_u*exp(exp_fact*sqrt(Uab/2*deltau));
    phi_d=phi_d*exp(-exp_fact*sqrt(Uab/2*deltau));
    Log_W=Log_W+0.5*shift^2-shift*exp_fact-sqrt(Uab*deltau/2)*n*exp_fact;
    for i=1:N+1
        if flag_bbp(i)==1
            B_u(list(i),i)=B_u(list(i),i)*exp(exp_fact*sqrt(Uab/2*deltau));
            B_d(list(i),i)=B_d(list(i),i)*exp(-exp_fact*sqrt(Uab/2*deltau));
        end
    end
end