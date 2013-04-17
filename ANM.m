% runs ANM code
clear
% close all

% type of diffusion coefficients
difftype = 'diffusion';
corr = 'no';

% load in xs data
reg1 = load('./extract_data/reg1.mat');
reg2 = load('./extract_data/reg2.mat');
% reg2.chi(:) = 0;
% reg1.keff = 1.02926;
% reg2.keff = 1.02926;

% load in form functions
reg1.form = load('./extract_data/reg1form.mat');
reg2.form = load('./extract_data/reg2form.mat');

% load in H1 correction curve
H1 = load('../correction_curve/PN_buckling/src/H-1600corr.mat');

% set lengths
reg1.L = 10.07872;
reg2.L = 10.07872;

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
% [reg1,reg2] = ANM_solve_fluxes(reg1,reg2);

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
dx1 = reg1.L/size(reg1.form.flux,2);
x1het = linspace(0 + dx1/2, reg1.L - dx1/2, size(reg1.form.flux,2));
dx2 = reg2.L/size(reg2.form.flux,2);
x2het = linspace(0 + dx2/2, reg2.L - dx2/2, size(reg2.form.flux,2));

figure(1)
%hold off
plot(x1,reg1.ANMphi1(x1));
hold on
plot(x2+reg1.L,reg2.ANMphi1(x2),'r');
title('Group 1 Flux');
plot(x1het,reg1.homflux(1,:),'k')
plot(x2het + reg1.L,reg2.homflux(1,:),'k') 


figure(2)
%hold off
plot(x1,reg1.ANMphi2(x1));
hold on
plot(x1+reg1.L,reg2.ANMphi2(x2),'r');
title('Group 2 Flux');
plot(x1het,reg1.homflux(2,:),'k')
plot(x2het + reg1.L,reg2.homflux(2,:),'k')

figure(3)
plot(reg1.MCpinpower)
hold on
plot(reg1.pinpower,'k.')
plot(9:16,reg2.pinpower,'k.')

figure(4)
plot(x1het,reg1.ANMphi1(x1het)./reg1.ANMphi2(x1het),'r')
hold on
plot(x2het+reg1.L,reg2.ANMphi1(x2het)./reg2.ANMphi2(x2het),'b')
plot([x1het,x2het + reg1.L],reg2.meshflux(1,:)./reg2.meshflux(2,:),'k')

