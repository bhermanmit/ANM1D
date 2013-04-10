function [reg1,reg2] = normalize_openmc(reg1,reg2)

% get x vectors
dx1 = reg1.L/size(reg1.form.flux,2);
x1 = linspace(0 + dx1/2, reg1.L - dx1/2, size(reg1.form.flux,2));
dx2 = reg2.L/size(reg2.form.flux,2);
x2 = linspace(0 + dx2/2, reg2.L - dx2/2, size(reg2.form.flux,2));

% integrate form functions
phireg1g1 = sum(reg1.form.flux(1,:));
phireg1g2 = sum(reg1.form.flux(2,:));
phireg2g1 = sum(reg2.form.flux(1,:));
phireg2g2 = sum(reg2.form.flux(2,:));

% adjust flux by integrated values
reg1.form.flux(1,:) = reg1.form.flux(1,:).*reg1.flux(1)./phireg1g1/dx1;
reg1.form.flux(2,:) = reg1.form.flux(2,:).*reg1.flux(2)./phireg1g2/dx1;
reg2.form.flux(1,:) = reg2.form.flux(1,:).*reg2.flux(1)./phireg2g1/dx2;
reg2.form.flux(2,:) = reg2.form.flux(2,:).*reg2.flux(2)./phireg2g2/dx2;

% integrate OpenMC distributions
phig1 = sum(reg1.meshflux(1,:));
phig2 = sum(reg1.meshflux(2,:));

% calculate normalization value
norm = (sum(reg1.flux) + sum(reg2.flux))/(phig1+phig2);

% adjust flux by normalization value
reg1.meshflux = reg1.meshflux*norm/dx1;
reg2.meshflux = reg2.meshflux*norm/dx1;

% integrate form functions
nsfreg1g1 = sum(reg1.form.nufission(1,:));
nsfreg1g2 = sum(reg1.form.nufission(2,:));
nsfreg2g1 = sum(reg2.form.nufission(1,:));
nsfreg2g2 = sum(reg2.form.nufission(2,:));

% adjust flux by integrated values
reg1.form.nufission(1,:) = reg1.form.nufission(1,:).*reg1.nsigf(1)*reg1.flux(1)./nsfreg1g1;
reg1.form.nufission(2,:) = reg1.form.nufission(2,:).*reg1.nsigf(2)*reg1.flux(2)./nsfreg1g2;
reg2.form.nufission(1,:) = reg2.form.nufission(1,:).*reg2.nsigf(1)*reg2.flux(1)./nsfreg2g1;
reg2.form.nufission(2,:) = reg2.form.nufission(2,:).*reg2.nsigf(2)*reg2.flux(2)./nsfreg2g2;

% integrate OpenMC distributions
nsfg1 = sum(reg1.meshnsf(1,:));
nsfg2 = sum(reg1.meshnsf(2,:));

% calculate normalization value
norm = (sum(reg1.nsigf.*reg1.flux) + sum(reg2.nsigf.*reg2.flux))/(nsfg1+nsfg2);

% adjust nufission by normalization value
reg1.meshnsf = reg1.meshnsf*norm;
reg2.meshnsf = reg2.meshnsf*norm;

end