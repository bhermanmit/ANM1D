% runs ANM code

% load in xs data
reg1 = load('./extract_data/reg1.mat');
reg2 = load('./extract_data/reg2.mat');

% set lengths
L1 = 10.3888;
L2 = 10.3888;

% set RDFs to 1
reg1.f(1:2) = 1;
reg2.f(1:2) = 1;

% compute neutron balance
[reg1,reg2] = neutron_balance(reg1,reg2);

% solve for fluxes
[reg1,reg2] = ANM_solve_fluxes(reg1,reg2,L1,L2);

% % plot fluxes
% x1 = linspace(0,L1,1000);
% x2 = linspace(0,L2,1000);
% figure(1)
% plot(x1,reg1.ANMphi1(x1));
% hold on
% plot(x1+L1,reg2.ANMphi1(x2),'r');
% title('Group 1 Flux');
% figure(2)
% plot(x1,reg1.ANMphi2(x1));
% hold on
% plot(x1+L1,reg2.ANMphi2(x2),'r');
% title('Group 2 Flux');

% compute discontinuity factors
[reg1,reg2] = ANM_compute_discontinuity(reg1,reg2,L1,L2);

% solve for fluxes
[reg1,reg2] = ANM_solve_fluxes(reg1,reg2,L1,L2);

% plot fluxes
x1 = linspace(0,L1,1000);
x2 = linspace(0,L2,1000);
figure(1)
plot(x1,reg1.ANMphi1(x1));
hold on
plot(x1+L1,reg2.ANMphi1(x2),'r');
title('Group 1 Flux');
figure(2)
plot(x1,reg1.ANMphi2(x1));
hold on
plot(x1+L1,reg2.ANMphi2(x2),'r');
title('Group 2 Flux');