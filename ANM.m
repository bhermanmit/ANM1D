% runs ANM code
clear
close all

% type of diffusion coefficients
difftype = 'transport';
corr = 'no';

% load in xs data
reg1 = load('./extract_data/reg1.mat');
reg2 = load('./extract_data/reg2.mat');

% load in form functions
reg1.form = load('./extract_data/reg1form.mat');
reg2.form = load('./extract_data/reg2form.mat');

% load in H1 correction curve
H1 = load('../correction_curve/PN_buckling/src/H-1corr.mat');

% set lengths
reg1.L = 10.3888;
reg2.L = 10.3888;

% set RDFs to 1
reg1.f(1:2) = 1;
reg2.f(1:2) = 1;

% compute neutron balance
[reg1,reg2] = neutron_balance(reg1,reg2);

% compute diffusion coefficients
[reg1,reg2] = compute_diffusion(reg1,reg2,H1,difftype,corr);

% normalize OpenMC distributions to average integrated fluxes
% [reg1,reg2] = normalize_openmc(reg1,reg2);

% solve for fluxes
[reg1,reg2] = ANM_solve_fluxes(reg1,reg2);

% compute discontinuity factors
[reg1,reg2] = ANM_compute_discontinuity(reg1,reg2);

% solve for fluxes
[reg1,reg2] = ANM_solve_fluxes(reg1,reg2);

% reconstruct fluxes
[reg1,reg2] = ANM_reconstruct_fluxes(reg1,reg2);

% reconstruct nu-fission
[reg1,reg2,err] = ANM_reconstruct_nufission(reg1,reg2);

% plot fluxes
x1 = linspace(0,reg1.L,1000);
x2 = linspace(0,reg2.L,1000);
figure(1)
plot(x1,reg1.ANMphi1(x1));
hold on
plot(x2+reg1.L,reg2.ANMphi1(x2),'r');
title('Group 1 Flux');


figure(2)
plot(x1,reg1.ANMphi2(x1));
hold on
plot(x1+reg1.L,reg2.ANMphi2(x2),'r');
title('Group 2 Flux');