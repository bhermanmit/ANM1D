function [reg1,reg2] = neutron_balance(reg1,reg2)

reg1.balance = reg1.curr + reg1.sigt.*reg1.flux - reg1.sigs'*reg1.flux - reg1.chi*reg1.nsigf'*reg1.flux/reg1.keff;
reg1.balance = reg1.balance./reg1.flux;

end