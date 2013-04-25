function [reg1,reg2] = compute_diffusion(reg1,reg2,H1,difftype,corr)

% last thermal group
cutoff = 47;

% compute fine distribution for transport rate
reg1.finetransrate = reg1.totalrate - reg1.p1rate;
reg2.finetransrate = reg2.totalrate - reg2.p1rate;

% correction fine distibution
if strcmp(corr,'yes')
    
    % calculate H-1 transport rate
    reg1.H1finetransrate = reg1.H1totalrate(1:46) - reg1.H1p1rate(1:46);
    reg2.H1finetransrate = reg2.H1totalrate(1:46) - reg2.H1p1rate(1:46);
    
    % remove this from fine transport rate
    reg1.finetransrate(1:46) = reg1.finetransrate(1:46) - reg1.H1finetransrate(1:46);
    reg2.finetransrate(1:46) = reg2.finetransrate(1:46) - reg2.H1finetransrate(1:46);
    
    % compute actual H1 fine transport rate
    reg1.H1finetransrate(1:46) = reg1.H1totalrate(1:46).*H1.H1rat(1:46);
    reg2.H1finetransrate(1:46) = reg2.H1totalrate(1:46).*H1.H1rat(1:46);
    
    % add H1 back to fine transport rate
    reg1.finetransrate(1:46) = reg1.finetransrate(1:46) + reg1.H1finetransrate(1:46);
    reg2.finetransrate(1:46) = reg2.finetransrate(1:46) + reg2.H1finetransrate(1:46);
    
end
% if strcmp(corr,'yes')
%     
%     % calculate H-1 transport rate
%     reg1.H1finetransrate = reg1.H1totalrate - reg1.H1p1rate;
%     reg2.H1finetransrate = reg2.H1totalrate - reg2.H1p1rate;
%     
%     % remove this from fine transport rate
%     reg1.finetransrate = reg1.finetransrate - reg1.H1finetransrate;
%     reg2.finetransrate = reg2.finetransrate - reg2.H1finetransrate;
%     
%     % compute actual H1 fine transport rate
%     reg1.H1finetransrate = reg1.H1totalrate.*H1.H1rat;
%     reg2.H1finetransrate = reg2.H1totalrate.*H1.H1rat;
%     
%     % add H1 back to fine transport rate
%     reg1.finetransrate = reg1.finetransrate + reg1.H1finetransrate;
%     reg2.finetransrate = reg2.finetransrate + reg2.H1finetransrate;
%     
% end

switch(difftype)
    
    case('transport')
        
        % calculate coarse distribution
        reg1.transrate(1,1) = sum(reg1.finetransrate(1:cutoff-1));
        reg1.transrate(2,1) = sum(reg1.finetransrate(cutoff:70));
        reg2.transrate(1,1) = sum(reg2.finetransrate(1:cutoff-1));
        reg2.transrate(2,1) = sum(reg2.finetransrate(cutoff:70));
        
        % calculate transport xs
        reg1.transxs = reg1.transrate./reg1.flux;
        reg2.transxs = reg2.transrate./reg2.flux;
        
        % overwrite diffusion coefficients
        reg1.diff = 1./(3*reg1.transxs);
        reg2.diff = 1./(3*reg2.transxs);
        
    case('diffusion')
        
        % calculate fine distribution of transport xs
        reg1.finetransxs = reg1.finetransrate./reg1.fluxrate;
        reg2.finetransxs = reg2.finetransrate./reg2.fluxrate;
        
        % calculate fine distribution of diffusion coeffs
        reg1.finediff = 1./(3*reg1.finetransxs);
        reg2.finediff = 1./(3*reg2.finetransxs);

        % calculate diffusion rate
        reg1.finediffrate = reg1.finediff.*reg1.fluxrate;
        reg2.finediffrate = reg2.finediff.*reg2.fluxrate;
        
        % collapse to coarse distribution
        reg1.diffrate(1,1) = sum(reg1.finediffrate(1:cutoff-1));
        reg1.diffrate(2,1) = sum(reg1.finediffrate(cutoff:70));
        reg2.diffrate(1,1) = sum(reg2.finediffrate(1:cutoff-1));
        reg2.diffrate(2,1) = sum(reg2.finediffrate(cutoff:70));
        
        % calculate transport xs
        reg1.diff = reg1.diffrate./reg1.flux;
        reg2.diff = reg2.diffrate./reg2.flux;
        
end

end