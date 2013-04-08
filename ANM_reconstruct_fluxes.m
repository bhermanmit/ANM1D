function [reg1,reg2] = ANM_reconstruct_fluxes(reg1,reg2)

% get x vectors
dx1 = reg1.L/size(reg1.form.flux,2);
x1 = linspace(0 + dx1/2, reg1.L - dx1/2, size(reg1.form.flux,2));
dx2 = reg2.L/size(reg2.form.flux,2);
x2 = linspace(0 + dx2/2, reg2.L - dx2/2, size(reg2.form.flux,2));

% synthesize
reg1.fluxrecon(1,:) = reg1.ANMphi1(x1).*reg1.form.flux(1,:);
reg1.fluxrecon(2,:) = reg1.ANMphi2(x1).*reg1.form.flux(2,:);
reg2.fluxrecon(1,:) = reg2.ANMphi1(x2).*reg2.form.flux(1,:);
reg2.fluxrecon(2,:) = reg2.ANMphi2(x2).*reg2.form.flux(2,:);

% renormalize to homogeneous integrated flux values
reg1.fluxrecon(1,:) = reg1.fluxrecon(1,:)*reg1.iphi1/sum(reg1.fluxrecon(1,:)*dx1);
reg1.fluxrecon(2,:) = reg1.fluxrecon(2,:)*reg1.iphi2/sum(reg1.fluxrecon(2,:)*dx1);
reg2.fluxrecon(1,:) = reg2.fluxrecon(1,:)*reg2.iphi1/sum(reg2.fluxrecon(1,:)*dx2);
reg2.fluxrecon(2,:) = reg2.fluxrecon(2,:)*reg2.iphi2/sum(reg2.fluxrecon(2,:)*dx2);

figure(5)
plot(x1,reg1.fluxrecon(1,:));
hold on
plot(x2 + reg1.L, reg2.fluxrecon(1,:), 'r');

figure(6)
plot(x1,reg1.fluxrecon(2,:));
hold on
plot(x2 + reg1.L, reg2.fluxrecon(2,:), 'r');


end