function [phi, O, w]=bos_norm(phi, Phi_T, N_wlk, O, w)
    new_O=zeros(N_wlk,1);
    new_w=zeros(N_wlk,1);
    for i=1:N_wlk
        new_phi=phi(:,i)/sqrt(phi(:,i)'*phi(:,i));
        new_O(i)=Phi_T'*new_phi;
        new_w(i)=w(i)*sqrt(phi(:,i)'*phi(:,i));
        phi(:,i)=new_phi;
    end
    O=new_O;
    w=new_w;
end