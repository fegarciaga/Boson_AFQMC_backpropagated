%% Initialize internal quantities
N_sites=Lx*Ly*Lz;
% form the one-body kinetic Hamiltonian
H_k=H_K(Lx, Ly, Lz, kx, ky, kz, tx, ty, tz);
% form the of-spin-sector matrix elements
J=jj*eye(N_sites);
% form joint-one-body operator 
K=[H_k-U/2*eye(N_sites) J; J H_k-U/2*eye(N_sites)]; 
%% Initialize the trial wave function and calculate the ensemble's initial energy 
% Diagonalize the one-body Hamiltonian to get the non-interacting single-particle orbitals:
[psi_nonint,E_nonint_m] = eig(K);
% assemble the non-interacting single-particle orbitals into a Permanent:
Phi_T=repmat(psi_nonint(:,1),1);

N_stps = 500;
 
% Performs Gross-Piaevskii mean field calculation
[Phi_T,E_T] = GP(Phi_T, N_stps, K, U, Uab, N_sites, N_par, deltau);

display(E_T);
display(Phi_T);

% Stores the mean field density profiles.
n_up = N_par*diag((Phi_T(1:N_sites,1)').'*Phi_T(1:N_sites,1).');
n_dn = N_par*diag((Phi_T(1+N_sites:2*N_sites,1)').'*Phi_T(1+N_sites:2*N_sites,1).');

% One body operator without mean field substraction is kept for measurements purposes
K_old=K;

% The one body operator is modified due to the mean field corrections
K(1:N_sites,1:N_sites)=K(1:N_sites,1:N_sites)+U*eye(N_sites,N_sites).*n_up+0.5*Uab*eye(N_sites,N_sites).*(n_up+n_dn)-0.5*Uab*eye(N_sites,N_sites).*(n_up-n_dn);
K(1+N_sites:2*N_sites,1+N_sites:2*N_sites)=K(1+N_sites:2*N_sites,1+N_sites:2*N_sites)+U*eye(N_sites,N_sites).*n_dn+0.5*Uab*eye(N_sites,N_sites).*(n_up+n_dn)+0.5*Uab*eye(N_sites,N_sites).*(n_up-n_dn);

% the matrix of the operator exp(-deltau*K/2)
Proj_k_half = expm(-0.5*deltau*K);

%% Assemble the initial population of walkers
Phi=zeros(2*N_sites,N_wlk);
% initiate each walker to be the trial wave function
for i=1:N_wlk
    % Phi(:,i) is the ith walker. Each is a matrix of size 2*N_sites by N_par
    % The left N_sites by N_up block is the spin up sector
    % The rest is the spin down sector
    % They are propagated independently and only share the auxiliary field
    Phi(:,i)=Phi_T; 
end

% initiate the weight and overlap of each walker to 1
w=ones(N_wlk,1);
O=ones(N_wlk,1);
% the arrays that store the energy and weight at each block
E_blk=zeros(N_blk,1);
W_blk=zeros(N_blk,1);
%% Prefactor due to one body elements
% fac_norm also include
fac_norm=real(E_T)*deltau+(0.5*U*real((n_up'*n_up)+(n_dn'*n_dn))+0.5*Uab*real(n_up'*n_dn+n_dn'*n_up))*deltau;
%% filename to be saved
savedFileName=strcat(int2str(Lx),'x',int2str(Ly),'x',int2str(Lz),'_',int2str(N_par),'par',num2str(U, '%4.2f'),'_kx',num2str(kx,'%+7.4f'),'_ky',num2str(ky,'%+7.4f'),'_kz',num2str(kz,'%+7.4f'),'_Nwlk_',int2str(N_wlk),suffix,'.mat');
%% randomize the random number generator seed based on the current time
rng('default');