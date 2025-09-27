function [alpha_hat,se_hat,sel,sel_rel,sel_irrel] = pdl_opt_lambda_sel_prob_tau(y,d,X,sig_y,sig_d,tau,s)
    
    [n,p] = size(X);

    yc = y - mean(y);
    Xc = X - mean(X);  
    dc = d - mean(d);
    
    lam_star_y = sig_y*2*sqrt(2*(1+tau)*log(p)/n);                
    beta_hat = lasso(Xc,yc,'Lambda',lam_star_y,'Standardize',true);
    lam_star_d = sig_d*2*sqrt(2*(1+tau)*log(p)/n);               
    gamma_hat = lasso(Xc,dc,'Lambda',lam_star_d,'Standardize',true);   
    Id = abs(gamma_hat)>0;
    Iy = abs(beta_hat)>0;
    I = max([Id,Iy],[],2);
    
    Xmat = [d X(:,I) ones(n,1)];
    b = (Xmat'*Xmat)\(Xmat'*y);
    VC = inv(Xmat'*Xmat)*sig_y.^2; 
    SE = sqrt(diag(VC));
    alpha_hat = b(1);
    se_hat = SE(1);
    
    sel = sum(I);
    
    Id_rel = abs(gamma_hat(1:s))>0;
    Iy_rel = abs(beta_hat(1:s))>0;
    I_rel = max([Id_rel,Iy_rel],[],2);
    sel_rel = sum(I_rel);
    
    Id_irrel = abs(gamma_hat((s+1):p))>0;
    Iy_irrel = abs(beta_hat((s+1):p))>0;
    I_irrel = max([Id_irrel,Iy_irrel],[],2);
    sel_irrel = sum(I_irrel);
    
end