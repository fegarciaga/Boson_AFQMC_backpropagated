function [phi] = bos_halfK(phi, Proj_k_half)
    %% propagate the walker by exp(-deltau*K/2)
    phi=Proj_k_half*phi;
end