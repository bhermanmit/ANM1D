% compute k-eff of two region two group problem

% load data
load mat

% Length of slabs
L1 = 10.3888;
L2 = 10.3888;

% compute removal cross sections
reg(1).sigr1 = mat(1).totxs(1) - mat(1).scattxs(1);
reg(1).sigr2 = mat(1).totxs(2) - mat(1).scattxs(4);
reg(2).sigr1 = mat(2).totxs(1) - mat(2).scattxs(1);
reg(2).sigr2 = mat(2).totxs(2) - mat(2).scattxs(4);

% extract downscatter and upscatter cross section
reg(1).sigs12 = mat(1).scattxs(2);
reg(2).sigs12 = mat(2).scattxs(2);
reg(1).sigs21 = mat(1).scattxs(3);
reg(2).sigs21 = mat(2).scattxs(3);

% extract nu-fission cross section
reg(1).nsigf1 = mat(1).nfissxs(1);
reg(1).nsigf2 = mat(1).nfissxs(3);
reg(2).nsigf1 = mat(2).nfissxs(1);
reg(2).nsigf2 = mat(2).nfissxs(3);

% extract diffusion coefficients
reg(1).D1 = mat(1).diff(1);
reg(1).D2 = mat(1).diff(2);
reg(2).D1 = mat(2).diff(1);
reg(2).D2 = mat(2).diff(2);

% set discontinuity factors
reg(1).f1 = 1.0189;
reg(1).f2 = 1.0398;
reg(2).f1 = 0.9811;
reg(2).f2 = 0.9602;
reg(1).f1 = 1.0;
reg(1).f2 = 1.0;
reg(2).f1 = 1.0;
reg(2).f2 = 1.0;

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

% plot fluxes
x1 = linspace(0,L1,1000);
x2 = linspace(0,L2,1000);
figure(3)
hold off
plot(x1,reg(1).phi1(x1));
hold on
plot(x1+L1,reg(2).phi1(x2),'r');
title('Group 1 Flux');
figure(4)
hold off
plot(x1,reg(1).phi2(x1));
hold on
plot(x1+L1,reg(2).phi2(x2),'r');
title('Group 2 Flux');