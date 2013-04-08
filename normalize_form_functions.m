function [reg1,reg2] = normalize_form_functions(reg1,reg2)

% get x vectors
dx1 = reg1.L/size(reg1.form.flux,2);
x1 = linspace(0 + dx1/2, reg1.L - dx1/2, size(reg1.form.flux,2));
dx2 = reg2.L/size(reg2.form.flux,2);
x2 = linspace(0 + dx2/2, reg2.L - dx2/2, size(reg2.form.flux,2));

% integrate form functions
phireg1g1 = trapz(x1,reg1.form.flux(1,:));
phireg1g2 = trapz(x1,reg1.form.flux(2,:));
phireg2g1 = trapz(x2,reg2.form.flux(1,:));
phireg2g2 = trapz(x2,reg2.form.flux(2,:));

% adjust flux by integrated values
reg1.form.flux(1,:) = reg1.form.flux(1,:).*reg1.flux(1)./phireg1g1;
reg1.form.flux(2,:) = reg1.form.flux(2,:).*reg1.flux(2)./phireg1g2;
reg2.form.flux(1,:) = reg2.form.flux(1,:).*reg2.flux(1)./phireg2g1;
reg2.form.flux(2,:) = reg2.form.flux(2,:).*reg2.flux(2)./phireg2g2;

figure(3)
plot(x1,reg1.form.flux(1,:));
hold on
plot(x2 + reg1.L, reg2.form.flux(1,:), 'r');

figure(4)
plot(x1,reg1.form.flux(2,:));
hold on
plot(x2 + reg1.L, reg2.form.flux(2,:), 'r');

end