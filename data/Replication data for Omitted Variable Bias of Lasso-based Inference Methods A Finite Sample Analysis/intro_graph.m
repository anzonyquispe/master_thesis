%% Preliminaries
% This code for reproducing the motivating figure in the introduction
% Authors: Kaspar Wuthrich and Ying Zhu
% DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
% Questions/error reports: kwuthrich@ucsd.edu

clear;
cd '/Users/kasparwuthrich/Dropbox/research/Kaspar_Ying_research/submission/ReStat/final submission/replication_package_final'
rng('default');

%% Setting

iter = 10000;
n = 500;
p = 200;
k = 5;

tau = 0.1;

c_vec = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];

%% Simulation

res_overall_bias = nan(length(c_vec),1);
res_overall_sel = nan(length(c_vec),1);

for j = 1:length(c_vec)

    res_bias = nan(iter,1);
    res_sel = nan(iter,1);

    c = c_vec(j);
    
    for it = 1:iter
        
        tbeta = zeros(p,1);
        tbeta(1:k) = 1;

        tgamma = zeros(p,1);
        tgamma(1:k) = 1;

        gamma = c*tgamma;
        beta = c*tbeta;

        X = randn(n,p);
        d = X*gamma + randn(n,1);
        y = X*beta + randn(n,1); % set alpha=0
            
        [alpha_hat,se_hat,sel,sel_rel,sel_irrel] = pdl_opt_lambda_sel_prob_tau(y,d,X,1,1,tau,k);
        
        res_bias(it,:) = alpha_hat;
        res_sel(it,:) = sel;
                
    end
        
    res_overall_bias(j,:) = mean(res_bias);
    res_overall_sel(j,:) = mean(res_sel);
    
end


%% Graphs

intro_illu = figure;
left_color = [0 0 0];
right_color = [0.5 0.5 0.5];
set(intro_illu,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left;
plot(c_vec,res_overall_bias,'-o','LineWidth',3);
ylabel('Bias');
ylim([0 .25]);

yyaxis right;
plot(c_vec,res_overall_sel,':o','LineWidth',3);
xlabel('Magnitude of \beta_j and \gamma_j for j=1,...,5');
ylim([0 6]);
xlim([min(c_vec) max(c_vec)]);
xticks(c_vec);
ylabel('Number of selected controls');

legend('Bias of post double Lasso','Number of selected controls','Location','northwest');

set(findall(gcf,'type','axes'),'fontsize',14);
set(findall(gcf,'type','text'),'fontSize',16);
set(findall(gcf,'type','legend'),'fontSize',16);
pbaspect([1.5 1 1]);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print('intro_illu','-dpdf','-bestfit');

close(intro_illu);


