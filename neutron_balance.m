function [reg1,reg2] = neutron_balance(reg1,reg2)

% compute efficientive downscatter
r1tmpsigs12 = (reg1.sigs(1,2)*reg1.flux(1) - reg1.sigs(2,1)*reg1.flux(2))/reg1.flux(1);
r2tmpsigs12 = (reg2.sigs(1,2)*reg2.flux(1) - reg2.sigs(2,1)*reg2.flux(2))/reg2.flux(1);
reg1.sigt(1) = reg1.sigt(1) - reg1.sigs(1,2) + r1tmpsigs12;
reg2.sigt(1) = reg2.sigt(1) - reg2.sigs(1,2) + r2tmpsigs12;
reg1.sigt(2) = reg1.sigt(2) - reg1.sigs(2,1);
reg2.sigt(2) = reg2.sigt(2) - reg2.sigs(2,1);
reg1.sigs(2,1) = 0;
reg2.sigs(2,1) = 0;
reg1.sigs(1,2) = r1tmpsigs12;
reg2.sigs(1,2) = r2tmpsigs12;


% region 1 balance
reg1.balance = (reg1.currR - reg1.currL) + reg1.sigt.*reg1.flux - reg1.sigs'*reg1.flux - reg1.chi*reg1.nsigf'*reg1.flux/reg1.keff;
reg1.balance = reg1.balance./reg1.flux;

% region 2 balance
reg2.balance = (reg2.currR - reg2.currL) + reg2.sigt.*reg2.flux - reg2.sigs'*reg2.flux - reg2.chi*reg2.nsigf'*reg2.flux/reg2.keff;
reg2.balance = reg2.balance./reg2.flux;

end