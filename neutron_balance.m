function [reg1,reg2] = neutron_balance(reg1,reg2)

reg1.balance = reg1.curr + reg1.sigt.*reg1.flux - reg1.sigs'*reg1.flux - reg1.chi*reg1.nsigf'*reg1.flux/reg1.keff;
reg1.balance = reg1.balance./reg1.flux;

reg2.balance = reg2.curr + reg2.sigt.*reg2.flux - reg2.sigs'*reg2.flux - reg2.chi*reg1.nsigf'*reg2.flux/reg2.keff;
reg2.balance = reg2.balance./reg2.flux;

end