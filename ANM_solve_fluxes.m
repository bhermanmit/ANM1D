function [reg1,reg2] = ANM_solve_fluxes(reg1,reg2)

% geometry
L1 = reg1.L;
L2 = reg2.L;

% compute removal cross sections
reg(1).sigr1 = reg1.sigt(1) - reg1.sigs(1,1);
reg(1).sigr2 = reg1.sigt(2) - reg1.sigs(2,2);
reg(2).sigr1 = reg2.sigt(1) - reg2.sigs(1,1);
reg(2).sigr2 = reg2.sigt(2) - reg2.sigs(2,2);

% set discontinuity factors
reg(1).f1 = reg1.f(1);
reg(1).f2 = reg1.f(2);
reg(2).f1 = reg2.f(1);
reg(2).f2 = reg2.f(2);

% resave other data in the structure below
reg(1).sigs12 = reg1.sigs(1,2);
reg(1).sigs21 = reg1.sigs(2,1);
reg(1).nsigf1 = reg1.nsigf(1);
reg(1).nsigf2 = reg1.nsigf(2);
reg(1).D1 = reg1.diff(1);
reg(1).D2 = reg1.diff(2);
reg(2).sigs12 = reg2.sigs(1,2);
reg(2).sigs21 = reg2.sigs(2,1);
reg(2).nsigf1 = reg2.nsigf(1);
reg(2).nsigf2 = reg2.nsigf(2);
reg(2).D1 = reg2.diff(1);
reg(2).D2 = reg2.diff(2);

% create lambda functions
reg(1).lamb1 = @(k) -reg(1).sigr2/(2*reg(1).D2) ...
    - reg(1).sigr1/(2*reg(1).D1) ...
    + reg(1).nsigf1/(2*reg(1).D1*k) ...
    + (reg(1).sigr2^2/(4*reg(1).D2^2) ...
    - reg(1).sigr1*reg(1).sigr2/(2*reg(1).D1*reg(1).D2) ...
    + reg(1).sigs12*reg(1).sigs21/(reg(1).D1*reg(1).D2) ...
    + reg(1).nsigf1*reg(1).sigr2/(2*reg(1).D1*reg(1).D2*k) ...
    + reg(1).nsigf2*reg(1).sigs12/(reg(1).D1*reg(1).D2*k) ...
    + reg(1).sigr1^2/(4*reg(1).D1^2) ...
    - reg(1).nsigf1*reg(1).sigr1/(2*reg(1).D1^2*k) ...
    + reg(1).nsigf1^2/(4*reg(1).D1^2*k^2))^(0.5);

reg(1).lamb2 = @(k) -reg(1).sigr2/(2*reg(1).D2) ...
    - reg(1).sigr1/(2*reg(1).D1) ...
    + reg(1).nsigf1/(2*reg(1).D1*k) ...
    - (reg(1).sigr2^2/(4*reg(1).D2^2) ...
    - reg(1).sigr1*reg(1).sigr2/(2*reg(1).D1*reg(1).D2) ...
    + reg(1).sigs12*reg(1).sigs21/(reg(1).D1*reg(1).D2) ...
    + reg(1).nsigf1*reg(1).sigr2/(2*reg(1).D1*reg(1).D2*k) ...
    + reg(1).nsigf2*reg(1).sigs12/(reg(1).D1*reg(1).D2*k) ...
    + reg(1).sigr1^2/(4*reg(1).D1^2) ...
    - reg(1).nsigf1*reg(1).sigr1/(2*reg(1).D1^2*k) ...
    + reg(1).nsigf1^2/(4*reg(1).D1^2*k^2))^(0.5);

reg(2).lamb1 = @(k) -reg(2).sigr2/(2*reg(2).D2) ...
    - reg(2).sigr1/(2*reg(2).D1) ...
    + reg(2).nsigf1/(2*reg(2).D1*k) ...
    + (reg(2).sigr2^2/(4*reg(2).D2^2) ...
    - reg(2).sigr1*reg(2).sigr2/(2*reg(2).D1*reg(2).D2) ...
    + reg(2).sigs12*reg(2).sigs21/(reg(2).D1*reg(2).D2) ...
    + reg(2).nsigf1*reg(2).sigr2/(2*reg(2).D1*reg(2).D2*k) ...
    + reg(2).nsigf2*reg(2).sigs12/(reg(2).D1*reg(2).D2*k) ...
    + reg(2).sigr1^2/(4*reg(2).D1^2) ...
    - reg(2).nsigf1*reg(2).sigr1/(2*reg(2).D1^2*k) ...
    + reg(2).nsigf1^2/(4*reg(2).D1^2*k^2))^(0.5);

reg(2).lamb2 = @(k) -reg(2).sigr2/(2*reg(2).D2) ...
    - reg(2).sigr1/(2*reg(2).D1) ...
    + reg(2).nsigf1/(2*reg(2).D1*k) ...
    - (reg(2).sigr2^2/(4*reg(2).D2^2) ...
    - reg(2).sigr1*reg(2).sigr2/(2*reg(2).D1*reg(2).D2) ...
    + reg(2).sigs12*reg(2).sigs21/(reg(2).D1*reg(2).D2) ...
    + reg(2).nsigf1*reg(2).sigr2/(2*reg(2).D1*reg(2).D2*k) ...
    + reg(2).nsigf2*reg(2).sigs12/(reg(2).D1*reg(2).D2*k) ...
    + reg(2).sigr1^2/(4*reg(2).D1^2) ...
    - reg(2).nsigf1*reg(2).sigr1/(2*reg(2).D1^2*k) ...
    + reg(2).nsigf1^2/(4*reg(2).D1^2*k^2))^(0.5);

% create eigenvector matrices
reg(1).v11 = @(k) reg(1).lamb1(k) + reg(1).sigr2/reg(1).D2;
reg(1).v12 = @(k) reg(1).lamb2(k) + reg(1).sigr2/reg(1).D2;
reg(1).v21 = @(k) reg(1).sigs12/reg(1).D2;
reg(1).v22 = @(k) reg(1).sigs12/reg(1).D2;
reg(2).v11 = @(k) reg(2).lamb1(k) + reg(2).sigr2/reg(2).D2;
reg(2).v12 = @(k) reg(2).lamb2(k) + reg(2).sigr2/reg(2).D2;
reg(2).v21 = @(k) reg(2).sigs12/reg(2).D2;
reg(2).v22 = @(k) reg(2).sigs12/reg(2).D2;

% create implicit matrix on k
M  = @(k) [ [ ...
    reg(1).f1*reg(1).v11(k)*cos(sqrt(reg(1).lamb1(k))*L1), ...      % (1,1)
    reg(1).f1*reg(1).v12(k)*cosh(sqrt(-reg(1).lamb2(k))*L1), ...    % (1,2)
   -reg(2).f1*reg(2).v11(k)*cos(sqrt(reg(2).lamb1(k))*L2), ...      % (1,3)
   -reg(2).f1*reg(2).v12(k)*cosh(sqrt(-reg(2).lamb2(k))*L2)]; ...   % (1,4)
   [reg(1).f2*reg(1).v21(k)*cos(sqrt(reg(1).lamb1(k))*L1), ...      % (2,1)
    reg(1).f2*reg(1).v22(k)*cosh(sqrt(-reg(1).lamb2(k))*L1), ...    % (2,2)
   -reg(2).f2*reg(2).v21(k)*cos(sqrt(reg(2).lamb1(k))*L2), ...      % (2,3)
   -reg(2).f2*reg(2).v22(k)*cosh(sqrt(-reg(2).lamb2(k))*L2)]; ...   % (2,4)
   [reg(1).D1*reg(1).v11(k)*sqrt(reg(1).lamb1(k))*sin(sqrt(reg(1).lamb1(k))*L1), ...    % (3,1)
   -reg(1).D1*reg(1).v12(k)*sqrt(-reg(1).lamb2(k))*sinh(sqrt(-reg(1).lamb2(k))*L1), ... % (3,2)
    reg(2).D1*reg(2).v11(k)*sqrt(reg(2).lamb1(k))*sin(sqrt(reg(2).lamb1(k))*L2), ...    % (3,3)
   -reg(2).D1*reg(2).v12(k)*sqrt(-reg(2).lamb2(k))*sinh(sqrt(-reg(2).lamb2(k))*L2)];... % (3,4)
   [reg(1).D2*reg(1).v21(k)*sqrt(reg(1).lamb1(k))*sin(sqrt(reg(1).lamb1(k))*L1), ...    % (4,1)
   -reg(1).D2*reg(1).v22(k)*sqrt(-reg(1).lamb2(k))*sinh(sqrt(-reg(1).lamb2(k))*L1), ... % (4,2)
    reg(2).D2*reg(2).v21(k)*sqrt(reg(2).lamb1(k))*sin(sqrt(reg(2).lamb1(k))*L2), ...    % (4,3)
   -reg(2).D2*reg(2).v22(k)*sqrt(-reg(2).lamb2(k))*sinh(sqrt(-reg(2).lamb2(k))*L2)];... % (4,4)
   ];

% create determinant function
f = @(k) det(M(k));

% create iteration
kL = 0.1;
kU = 1.9;
kold = 1.0;
kLbase = kL;
kUbase = kU;
kountL = 0;
kountU = 0;
for i = 1:10000
    
    % evaluate function bounds
    fkU = f(kU);
    fkL = f(kL);
    
    % find root with secant method
    keff = kU - (fkU*(kU - kL))/(fkU - fkL);
    
    % evalue function of root
    fk = f(keff);
    
    % decide how to move brackets
    if fk*fkL > 0
        kL = keff;
    else
        kU = keff;
    end
    
    % check convergence errors
    err = abs(keff - kold);
    fprintf('ITER: %d  ERR: %d\n',i,err);
    if err < 1.e-8; break; end;
    kold = keff;
    
    % check to see if bound(s) are stagnant
    if (abs(kL - kLbase) < 1e-8)
        kountL = kountL + 1;
    end
    if (abs(kU - kUbase) < 1e-8)
        kountU = kountU + 1;
    end
    if (kountL > 2)
        kL = (kLbase - keff)/2;
        kLbase = kL;
        kountL = 0;
    end
    if (kountU > 2)
        kU = (kUbase - k)/2;
        kUbase = kU;
        kountU = 0;
    end
end % Modified False Position Method

% create matrix to solve for coeffs
A  = @(k) [ [ ...
    reg(1).f1*reg(1).v11(k)*cos(sqrt(reg(1).lamb1(k))*L1), ...      % (1,1)
    reg(1).f1*reg(1).v12(k)*cosh(sqrt(-reg(1).lamb2(k))*L1), ...    % (1,2)
   -reg(2).f1*reg(2).v11(k)*cos(sqrt(reg(2).lamb1(k))*L2)]; ...     % (1,3)
   [reg(1).f2*reg(1).v21(k)*cos(sqrt(reg(1).lamb1(k))*L1), ...      % (2,1)
    reg(1).f2*reg(1).v22(k)*cosh(sqrt(-reg(1).lamb2(k))*L1), ...    % (2,2)
   -reg(2).f2*reg(2).v21(k)*cos(sqrt(reg(2).lamb1(k))*L2)]; ...     % (2,3)
   [reg(1).D1*reg(1).v11(k)*sqrt(reg(1).lamb1(k))*sin(sqrt(reg(1).lamb1(k))*L1), ...    % (3,1)
   -reg(1).D1*reg(1).v12(k)*sqrt(-reg(1).lamb2(k))*sinh(sqrt(-reg(1).lamb2(k))*L1), ... % (3,2)
    reg(2).D1*reg(2).v11(k)*sqrt(reg(2).lamb1(k))*sin(sqrt(reg(2).lamb1(k))*L2)]; ...   % (3,3)
   ];

b = @(k) [ ...
    reg(2).f1*reg(2).v12(k)*cosh(sqrt(-reg(2).lamb2(k))*L2);
    reg(2).f2*reg(2).v22(k)*cosh(sqrt(-reg(2).lamb2(k))*L2);
    reg(2).D1*reg(2).v12(k)*sqrt(-reg(2).lamb2(k))*sinh(sqrt(-reg(2).lamb2(k))*L2)];

coeffs = A(keff)\b(keff);
reg(1).a = coeffs(1);
reg(1).c = coeffs(2);
reg(2).a = coeffs(3);
reg(2).c = 1;

% set up fluxes
reg(1).phi1 = @(x) reg(1).v11(keff)*reg(1).a*cos(sqrt(reg(1).lamb1(keff))*x) + ...
    reg(1).v12(keff)*reg(1).c*cosh(sqrt(-reg(1).lamb2(keff))*x);
reg(1).phi2 = @(x) reg(1).v21(keff)*reg(1).a*cos(sqrt(reg(1).lamb1(keff))*x) + ...
    reg(1).v22(keff)*reg(1).c*cosh(sqrt(-reg(1).lamb2(keff))*x);
reg(2).phi1 = @(x) reg(2).v11(keff)*reg(2).a*cos(sqrt(reg(2).lamb1(keff))*(x-L2)) + ...
    reg(2).v12(keff)*reg(2).c*cosh(sqrt(-reg(2).lamb2(keff))*(x-L2));
reg(2).phi2 = @(x) reg(2).v21(keff)*reg(2).a*cos(sqrt(reg(2).lamb1(keff))*(x-L2)) + ...
    reg(2).v22(keff)*reg(2).c*cosh(sqrt(-reg(2).lamb2(keff))*(x-L2));

% integral of fluxes
% reg(1).iphi1 = reg(1).v11(keff)*reg(1).a/sqrt(reg(1).lamb1(keff))*sin(sqrt(reg(1).lamb1(keff))*L1) + ...
%     reg(1).v12(keff)*reg(1).c/sqrt(-reg(1).lamb2(keff))*sinh(sqrt(-reg(1).lamb2(keff))*L1);
% reg(1).iphi2 = reg(1).v21(keff)*reg(1).a/sqrt(reg(1).lamb1(keff))*sin(sqrt(reg(1).lamb1(keff))*L1) + ...
%     reg(1).v22(keff)*reg(1).c/sqrt(-reg(1).lamb2(keff))*sinh(sqrt(-reg(1).lamb2(keff))*L1);
% reg(2).iphi1 = reg(2).v11(keff)*reg(2).a/sqrt(reg(2).lamb1(keff))*sin(sqrt(reg(2).lamb1(keff))*L2) + ...
%     reg(2).v12(keff)*reg(2).c/sqrt(-reg(2).lamb2(keff))*sinh(sqrt(-reg(2).lamb2(keff))*L2);
% reg(2).iphi2 = reg(2).v21(keff)*reg(2).a/sqrt(reg(2).lamb1(keff))*sin(sqrt(reg(2).lamb1(keff))*L2) + ...
%     reg(2).v22(keff)*reg(2).c/sqrt(-reg(2).lamb2(keff))*sinh(sqrt(-reg(2).lamb2(keff))*L2);

% integral of fluxes
reg(1).iphi1 = integral(reg(1).phi1,0,L1,'AbsTol',1e-16);
reg(1).iphi2 = integral(reg(1).phi2,0,L1,'AbsTol',1e-16);
reg(2).iphi1 = integral(reg(2).phi1,0,L2,'AbsTol',1e-16);
reg(2).iphi2 = integral(reg(2).phi2,0,L2,'AbsTol',1e-16);

% calculate normalization factor
reg1.flux = sum(reg2.meshflux(:,1:320),2);
reg2.flux = sum(reg2.meshflux(:,321:640),2);
norm = (sum(reg1.flux)+sum(reg2.flux))/(reg(1).iphi1+reg(1).iphi2+reg(2).iphi1+reg(2).iphi2);

% adjust fluxes
reg(1).phi1 = @(x) norm*reg(1).phi1(x);
reg(1).phi2 = @(x) norm*reg(1).phi2(x);
reg(2).phi1 = @(x) norm*reg(2).phi1(x);
reg(2).phi2 = @(x) norm*reg(2).phi2(x);
% reg(1).phi1 = @(x) reg(1).phi1(x)*reg1.flux(1)/reg(1).iphi1;
% reg(1).phi2 = @(x) reg(1).phi2(x)*reg1.flux(2)/reg(1).iphi2;
% reg(2).phi1 = @(x) reg(2).phi1(x)*reg2.flux(1)/reg(2).iphi1;
% reg(2).phi2 = @(x) reg(2).phi2(x)*reg2.flux(2)/reg(2).iphi2;

% save fluxes in output object
reg1.ANMphi1 = reg(1).phi1;
reg1.ANMphi2 = reg(1).phi2;
reg2.ANMphi1 = reg(2).phi1;
reg2.ANMphi2 = reg(2).phi2;
reg1.iphi1 = reg(1).iphi1*norm;
reg1.iphi2 = reg(1).iphi2*norm;
reg2.iphi1 = reg(2).iphi1*norm;
reg2.iphi2 = reg(2).iphi2*norm;
reg1.ANMkeff = keff;
reg2.ANMkeff = keff;

end
