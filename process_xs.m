function [reg1,reg2] = process_xs(reg1,reg2)

% compute removal xs
reg1.sigr(1) = reg1.sigt(1) - reg1.sigs(1,1);
reg1.sigr(2) = reg1.sigt(2) - reg1.sigs(2,2);
reg2.sigr(1) = reg2.sigt(1) - reg2.sigs(1,1);
reg2.sigr(2) = reg2.sigt(2) - reg2.sigs(2,2);

end