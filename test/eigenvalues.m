% symbolics
syms k
syms D1
syms D2
syms nsigf1
syms nsigf2
syms sigr1
syms sigr2
syms sigs12

% create matrix
F = [(nsigf1/k-sigr1)/D1, (nsigf2/k)/D1; sigs12/D2, -sigr2/D2];

% get eigenvalues and eigenvectors
[v,e] = eig(F);