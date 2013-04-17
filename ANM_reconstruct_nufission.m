function [reg1,reg2,err] = ANM_reconstruct_nufission(reg1,reg2)

% get x vectors
dx1 = reg1.L/8;
dx2 = reg2.L/8;

% compute nu-fission mesh
% for i = 1:8
%     reg1.nsfpinrate(:,i) = sum(reg1.form.nufission(:,1+(i-1)*20:i*20),2);
%     reg1.fluxpinrate(:,i) = sum(reg1.form.flux(:,1+(i-1)*20:i*20),2);
%     reg2.nsfpinrate(:,i) = sum(reg2.form.nufission(:,1+(i-1)*20:i*20),2);
%     reg2.fluxpinrate(:,i) = sum(reg2.form.flux(:,1+(i-1)*20:i*20),2);
%     
%     reg1.nsfMCpinrate(:,i) = sum(reg1.meshnsf(:,1+(i-1)*20:i*20),2);
%     j = i + 8;
%     reg1.nsfMCpinrate(:,j) = sum(reg1.meshnsf(:,1+(j-1)*20:j*20),2);
% end
% reg2.nsfMCpinrate = reg1.nsfMCpinrate;

for i = 1:8
    reg1.nsfpinrate(:,i) = sum(reg1.meshnsf(:,1+(i-1)*20:i*20),2);
    reg1.fluxpinrate(:,i) = sum(reg1.meshflux(:,1+(i-1)*20:i*20),2);
    j = i+8;
    reg2.nsfpinrate(:,i) = sum(reg2.meshnsf(:,1+(j-1)*20:j*20),2);
    reg2.fluxpinrate(:,i) = sum(reg2.meshflux(:,1+(j-1)*20:j*20),2);
    
    reg1.nsfMCpinrate(:,i) = sum(reg1.meshnsf(:,1+(i-1)*20:i*20),2);
    j = i + 8;
    reg1.nsfMCpinrate(:,j) = sum(reg1.meshnsf(:,1+(j-1)*20:j*20),2);
end
reg2.nsfMCpinrate = reg1.nsfMCpinrate;

% compute integral
reg1.nsf = sum(reg1.nsfMCpinrate(:,1:8),2);
reg2.nsf = sum(reg2.nsfMCpinrate(:,9:16),2);

% compute xs
reg1.nsfxs = reg1.nsfpinrate./reg1.fluxpinrate;
reg2.nsfxs = reg2.nsfpinrate./reg2.fluxpinrate;

% compute pin flux rates
for i = 1:8
    reg1.hompinflux(1,i) = quad(reg1.ANMphi1,(i-1)*dx1,i*dx1,1e-16);
    reg1.hompinflux(2,i) = quad(reg1.ANMphi2,(i-1)*dx1,i*dx1,1e-16);
    reg2.hompinflux(1,i) = quad(reg2.ANMphi1,(i-1)*dx2,i*dx2,1e-16);
    reg2.hompinflux(2,i) = quad(reg2.ANMphi2,(i-1)*dx2,i*dx2,1e-16);
end

% compute reconstructed nsf
reg1.nsfpinrecon = reg1.hompinflux.*reg1.nsfxs;
reg2.nsfpinrecon = reg2.hompinflux.*reg2.nsfxs;

% renormalize
reg1.nsfpinrecon(1,:) = reg1.nsfpinrecon(1,:)*reg1.nsf(1)/sum(reg1.nsfpinrecon(1,:));
reg1.nsfpinrecon(2,:) = reg1.nsfpinrecon(2,:)*reg1.nsf(2)/sum(reg1.nsfpinrecon(2,:));
reg2.nsfpinrecon(1,:) = reg2.nsfpinrecon(1,:)*reg2.nsf(1)/sum(reg2.nsfpinrecon(1,:));
reg2.nsfpinrecon(2,:) = reg2.nsfpinrecon(2,:)*reg2.nsf(2)/sum(reg2.nsfpinrecon(2,:));

% get pin powers
reg1.pinpower = sum(reg1.nsfpinrecon,1);
reg2.pinpower = sum(reg2.nsfpinrecon,1);

% get MC pin powers
reg1.MCpinpower = sum(reg1.nsfMCpinrate,1);
reg2.MCpinpower = sum(reg2.nsfMCpinrate,1);

% renormalize powers about 1
norm = 16/sum(reg1.pinpower + reg2.pinpower);
reg1.pinpower = reg1.pinpower*norm;
reg2.pinpower = reg2.pinpower*norm;
reg1.MCpinpower = reg1.MCpinpower*16/sum(reg1.MCpinpower);
reg2.MCpinpower = reg2.MCpinpower*16/sum(reg2.MCpinpower);

% computer L2 norm error
diff = [reg1.pinpower,reg2.pinpower] - reg1.MCpinpower;
err = sqrt((1/16)*sum(diff.^2))*100;

end