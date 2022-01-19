%  A script to set the input parameters and run a CPMC calculation
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% �2014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% system parameters:
Lx=3; % The number of lattice sites in the x direction
Ly=1; % The number of lattice sites in the y direction
Lz=1; % The number of lattice sites in the z direction

N_par=3;% The number of bosons

kx=0.0; % The x component of the twist angle in TABC (twist-averaging boundary condition)
ky=+0.0; % The y component of the twist angle in TABC
kz=0; % The z component of the twist angle in TABC

U=0; % The on-site repulsion strength in the Hubbard Hamiltonian
Uab=[8 10 12];
tx=1; % The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; % The hopping amplitude between nearest neighbor sites in the y direction
tz=1; % The hopping amplitude between nearest neighbor sites in the z direction
jj=0.1;

%% run parameters:
deltau=0.01; % The imaginary time step
N_wlk=200; % The number of random walkers
N_blksteps=40; % The number of random walk steps in each block
N_eqblk=50; % The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=90; % The number of blocks used in the measurement phase
itv_pc=5; % The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_nrm=1;
itv_Em=40; % The interval between two adjacent energy measurements
t_bp=121;
t_pop=1;
suffix=datestr(now,'_yymmdd_HHMMSS'); % time stamp for the saved *.mat filename. Can be changed to any desired string 

N_run1=length(Uab);
N_run2=3;

E_ave=zeros(N_run1,1);
E_err=zeros(N_run1,1);
D_ave=zeros(N_run1,2*Lx);
D_err=zeros(N_run1,2*Lx);
D_ave1=zeros(N_run1,2*Lx);
D_err1=zeros(N_run1,2*Lx);

for ii=1:N_run1
    suffix=strcat('_try',int2str(Uab(ii)));
    [E_ave(ii),E_err(ii),D_ave(ii,:),D_err(ii,:),savedFile]=PPMC_Bos(Lx,Ly,Lz,Lx,kx,ky,kz,U,Uab(ii),jj,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_pc,itv_nrm,itv_Em,t_bp,t_pop,suffix);
    [E_ave(ii),E_err(ii),D_ave1(ii,:),D_err1(ii,:),savedFile]=PPMC_Bos_w_prof(Lx,Ly,Lz,Lx,kx,ky,kz,U,Uab(ii),jj,tx,ty,tz,deltau,D_ave(ii,:),N_wlk,N_blksteps,N_eqblk,N_blk,itv_pc,itv_nrm,itv_Em,t_bp,t_pop,suffix);
end

%% post-run:
% load saved data into workspace for post-run analysis:
%load (savedFile);
%% Explanation of saved quantities:
% E: the array of energy of each block
% time: The total computational time
% E_nonint_v: the non-interacting energy levels of the system
% Phi_T: the trial wave function
% For other saved quantities, type "help CPMC_Lab"