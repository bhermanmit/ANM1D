% runs ANM code
clear
% close all

% type of diffusion coefficients
difftype = 'diffusion';
corr = 'yes';

% load in xs data
reg1 = load('./extract_data/reg1.mat');
reg2 = load('./extract_data/reg2.mat');
reg1SA = load('./extract_data/reg1SA.mat');
reg2SA = load('./extract_data/reg2SA.mat');
reg2.chi(:) = 0;
reg2SA.chi(:) = 0;
reg1.keff = 1.03035;
reg2.keff = 1.03035;
reg1SA.keff = 1.03035;
reg2SA.keff = 1.03035;

% load in H1 correction curve
H1 = load('../correction_curve/PN_buckling/src/H-1600corr.mat');

% set lengths
reg1.L = 10.07872*2;
reg2.L = 10.07872*2;
reg1SA.L = 10.07872*2;
reg2SA.L = 10.07872*2;

% set RDFs to 1
reg1.f(1:2) = 1;
reg2.f(1:2) = 1;

% compute neutron balance
[reg1,reg2] = neutron_balance(reg1,reg2);
[reg1SA,reg2SA] = neutron_balance(reg1SA,reg2SA);

% compute diffusion coefficients
[reg1,reg2] = compute_diffusion(reg1,reg2,H1,difftype,corr);
[reg1SA,reg2SA] = compute_diffusion(reg1SA,reg2SA,H1,difftype,corr);

% compute discontinuity factors
[reg1,reg2] = ANM_compute_discontinuity(reg1,reg2,'adfs');
% reg1.f = [9.905986239275264e-01     1.048814831170608e+00];
% reg2.f = [1 1];
reg1.f = [1 1];
reg2.f = [1.110858237865663e+00     8.385291258745322e-01];

% compute reference solution
[reg1,reg2] = ANM_solve_fluxes(reg1,reg2);

% copy over discontinuity factors
[reg1,reg2] = ANM_compute_discontinuity(reg1,reg2,'adfs');
reg2SA.f = reg2.f;
%reg1SA.f = reg1.f;
[reg1SA,reg2] = ANM_compute_discontinuity(reg1SA,reg2,'adfs');
%reg2SA.f = reg2.f;
reg1SA.f = [1 1];
%reg2SA.f = [1.110858237865663e+00     8.385291258745322e-01];
%reg2SA.f = [1 1];

% solve for fluxes
[regSA1,reg2SA] = ANM_solve_fluxes(reg1SA,reg2SA);

% plot fluxes
x1 = linspace(0,reg1.L,1000);
x2 = linspace(0,reg2.L,1000);
dx1 = reg1.L/(size(reg1.meshflux,2)/2);
x1het = linspace(0 + dx1/2, reg1.L - dx1/2, size(reg2.meshflux,2)/2);
dx2 = reg2.L/(size(reg2.meshflux,2)/2);
x2het = linspace(0 + dx2/2, reg2.L - dx2/2, size(reg2.meshflux,2)/2);

figure(1)
%hold off
plot(x1,reg1SA.ANMphi1(x1),'b--');
hold on
plot(x2+reg1.L,reg2SA.ANMphi1(x2),'r--');
title('Group 1 Flux');
plot([x1het,x2het + reg1.L],reg2.meshflux(1,:)/dx1,'k')

figure(2)
%hold off
plot(x1,reg1SA.ANMphi2(x1),'b--');
hold on
plot(x1+reg1.L,reg2SA.ANMphi2(x2),'r--');
title('Group 2 Flux');
plot([x1het,x2het + reg1.L],reg2.meshflux(2,:)/dx2,'k')