%% Preliminaries
% Code for reproducing the numerical example on the performance of post double Lasso
% Authors: Kaspar Wuthrich and Ying Zhu
% DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
% Questions/error reports: kwuthrich@ucsd.edu

clear;
cd '/Users/kasparwuthrich/Dropbox/research/Kaspar_Ying_research/submission/ReStat/final submission/replication_package_final'

rng('default');

%% Simulations

% Setup
iter = 50000;
n = 500;
p = 200;
k = 5;

Rsq1 = 0.5;
Rsq2 = 0.3;
Rsq3 = 0.1;
Rsq4 = 0.01;

tau = 0.1;

% Simulations

res_alpha = nan(iter,4);
res_sel = nan(iter,4);
res_sel_irrel = nan(iter,4);

for it = 1:iter
        
    [y1,d1,X1] = gen_design_r2(n,p,k,Rsq1,Rsq1);
    [y2,d2,X2] = gen_design_r2(n,p,k,Rsq2,Rsq2);
    [y3,d3,X3] = gen_design_r2(n,p,k,Rsq3,Rsq3);
    [y4,d4,X4] = gen_design_r2(n,p,k,Rsq4,Rsq4);

    [alpha1,se1,sel1,sel_rel1,sel_irrel1] = pdl_opt_lambda_sel_prob_tau(y1,d1,X1,1,1,tau,k);
    [alpha2,se2,sel2,sel_rel2,sel_irrel2] = pdl_opt_lambda_sel_prob_tau(y2,d2,X2,1,1,tau,k);
    [alpha3,se3,sel3,sel_rel3,sel_irrel3] = pdl_opt_lambda_sel_prob_tau(y3,d3,X3,1,1,tau,k);
    [alpha4,se4,sel4,sel_rel4,sel_irrel4] = pdl_opt_lambda_sel_prob_tau(y4,d4,X4,1,1,tau,k);

    res_alpha(it,:) = [alpha1 alpha2 alpha3 alpha4];
    res_sel(it,:) = [sel_rel1 sel_rel2 sel_rel3 sel_rel4];
    res_sel_irrel(it,:) = [sel_irrel1 sel_irrel2 sel_irrel3 sel_irrel4];

end

% Check that the selection probability is small for irrelevant regressors

mean(res_sel_irrel)

% Check whether distributions are skewed

kurtosis(res_alpha)

%% Figures

xmin = -.2; % lower bound
xmax = 0.4; % upper bound
h = 0.02; % bin width
edges = xmin:h:xmax;

fsd_dml = figure;
set(gcf, 'color','None');
set(gca, 'color','None');

subplot(2,2,1);
histogram(res_alpha(:,1),edges,'FaceColor',[0.5 0.5 0.5]);
hold;
f = @(z) iter*h*normpdf(z,0,1/sqrt(500));
fplot(f,[-0.2,0.4],'color',[0 0 0],'LineWidth',1.5); 
hold off;
xline(0,'--','True value, \alpha^*','LineWidth',1.5); 
xlim([-0.2,0.4]);
ylim([0 10000]); 
ylabel(''); 
title('(a) R^2 = 0.5');
set(gca,'YTick',[]);

subplot(2,2,2);
histogram(res_alpha(:,2),edges,'FaceColor',[0.5 0.5 0.5]);
hold;
f = @(z) iter*h*normpdf(z,0,1/sqrt(500));
fplot(f,[-0.2,0.4],'color',[0 0 0],'LineWidth',1.5); 
hold off;
xline(0,'--','True value, \alpha^*','LineWidth',1.5); 
xlim([-0.2,0.4]);
ylim([0 10000]); 
ylabel('');
title('(b) R^2 = 0.3');
set(gca,'YTick',[]);

subplot(2,2,3);
histogram(res_alpha(:,3),edges,'FaceColor',[0.5 0.5 0.5]);
hold;
f = @(z) iter*h*normpdf(z,0,1/sqrt(500));
fplot(f,[-0.2,0.4],'color',[0 0 0],'LineWidth',1.5); 
hold off;
xline(0,'--','True value, \alpha^*','LineWidth',1.5); 
xlim([-0.2,0.4]);
ylim([0 10000]); 
ylabel(''); 
title('(c) R^2 = 0.1');
set(gca,'YTick',[]);

subplot(2,2,4);
histogram(res_alpha(:,4),edges,'FaceColor',[0.5 0.5 0.5]);
hold;
f = @(z) iter*h*normpdf(z,0,1/sqrt(500));
fplot(f,[-0.2,0.4],'color',[0 0 0],'LineWidth',1.5); 
hold off;
xline(0,'--','True value, \alpha^*','LineWidth',1.5); 
xlim([-0.2,0.4]);
ylim([0 10000]); 
ylabel(''); 
title('(d) R^2 = 0.01');
set(gca,'YTick',[]);

set(findall(gcf,'type','axes'),'fontsize',10);
set(findall(gcf,'type','text'),'fontSize',10);
set(findall(gcf,'type','legend'),'fontSize',10);

axesHandles = findobj(get(fsd_dml,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print('fsd_dml','-dpdf','-bestfit');
close(fsd_dml);

xmin = -0.5; % lower bound
xmax = 5.5; % upper bound
h = 1; % bin width
edges = xmin:h:xmax;

hist_sel_dml = figure;
set(gcf, 'color','None');
set(gca, 'color','None');

subplot(2,2,1);
histogram(res_sel(:,1),edges,'FaceColor',[0.5 0.5 0.5]);
xlim([-0.5 5.5]);
ylim([0 50000]); 
title('(a) R^2 = 0.5');
xticks([0 1 2 3 4 5]);
ylim([0 50000]); 
set(gca,'YTick',[]);

subplot(2,2,2);
histogram(res_sel(:,2),edges,'FaceColor',[0.5 0.5 0.5]);
xlim([-0.5 5.5]);
ylim([0 50000]); 
ylabel(''); 
title('(b) R^2 = 0.3');
xticks([0 1 2 3 4 5]);
set(gca,'YTick',[]);

subplot(2,2,3);
histogram(res_sel(:,3),edges,'FaceColor',[0.5 0.5 0.5]);
xlim([-0.5 5.5]);
ylim([0 50000]); 
title('(c) R^2 = 0.1');
xticks([0 1 2 3 4 5]);
set(gca,'YTick',[]);

subplot(2,2,4);
histogram(res_sel(:,4),edges,'FaceColor',[0.5 0.5 0.5]);
xlim([-0.5 5.5]);
ylim([0 50000]); 
title('(d) R^2 = 0.01');
xticks([0 1 2 3 4 5]);
set(gca,'YTick',[]);

set(findall(gcf,'type','axes'),'fontsize',10);
set(findall(gcf,'type','text'),'fontSize',10);
set(findall(gcf,'type','legend'),'fontSize',10);

axesHandles = findobj(get(hist_sel_dml,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print('hist_sel_dml','-dpdf','-bestfit');
close(hist_sel_dml);


%% Simulation large sample

% Setup
iter = 2000;
n = 14238;
p = 384;
k = 10;
tau = 0.1;

Rsq1 = 0.2;
Rsq2 = 0.01;

% Simulations

res_alpha_large = nan(iter,2);
res_sel_large = nan(iter,2);
res_sel_irrel_large = nan(iter,2);

for it = 1:iter
        
    [y1,d1,X1] = gen_design_r2(n,p,k,Rsq1,Rsq1);
    [y2,d2,X2] = gen_design_r2(n,p,k,Rsq2,Rsq2);

    [alpha1,se1,sel1,sel_rel1,sel_irrel1] = pdl_opt_lambda_sel_prob_tau(y1,d1,X1,1,1,tau,k);
    [alpha2,se2,sel2,sel_rel2,sel_irrel2] = pdl_opt_lambda_sel_prob_tau(y2,d2,X2,1,1,tau,k);

    res_alpha_large(it,:) = [alpha1 alpha2];
    res_sel_large(it,:) = [sel_rel1 sel_rel2];
    res_sel_irrel_large(it,:) = [sel_irrel1 sel_irrel2];

end


% Check that the selection probability is small for irrelevant regressors

mean(res_sel_irrel_large)

%% Figures

xmin = -0.0525; % lower bound
xmax = 0.1; % upper bound
h = 0.005; % bin width
edges = xmin:h:xmax;

xmin2 = -0.5; % lower bound
xmax2 = 10.5; % upper bound
h2 = 1; % bin width
edges2 = xmin2:h2:xmax2;

fsd_large = figure('Position',[1000,1000,625,251]);
set(gcf, 'color','None');
set(gca, 'color','None');

subplot(1,2,1);
histogram(res_alpha_large(:,1),edges,'FaceColor',[0.5 0.5 0.5]);
hold;
f = @(z) iter*h*normpdf(z,0,1/sqrt(n));
fplot(f,[-0.05 0.075],'color',[0 0 0],'LineWidth',1.5); 
hold off
xline(0,'--','True value, \alpha^*','LineWidth',1.5); 
xlim([-0.05 0.075]);
ylim([0 550]); 
ylabel(''); 
title('R^2 = 0.2');
set(gca,'YTick',[]);

subplot(1,2,2);
histogram(res_alpha_large(:,2),edges,'FaceColor',[0.5 0.5 0.5]);
hold;
f = @(z) iter*h*normpdf(z,0,1/sqrt(n));
fplot(f,[-0.05 0.1],'color',[0 0 0],'LineWidth',1.5); 
hold off;
xline(0,'--','True value, \alpha^*','LineWidth',1.5); 
xlim([-0.05 0.075]);
ylim([0 550]); 
ylabel(''); 
title('R^2 = 0.01');
set(gca,'YTick',[]);

set(findall(gcf,'type','axes'),'fontsize',12);
set(findall(gcf,'type','text'),'fontSize',12);
set(findall(gcf,'type','legend'),'fontSize',12);
axesHandles = findobj(get(fsd_large,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print('fsd_large','-dpdf','-bestfit');
close(fsd_large);

hist_sel_large = figure('Position',[1000,1000,625,251]);
set(gcf, 'color','None');
set(gca, 'color','None');

subplot(1,2,1);
histogram(res_sel_large(:,1),edges2,'FaceColor',[0.5 0.5 0.5]);
xlim([-0.5 10.5]);
ylim([0 2000]); 
title('R^2 = 0.2');
set(gca,'YTick',[]);

subplot(1,2,2);
histogram(res_sel_large(:,2),edges2,'FaceColor',[0.5 0.5 0.5]);
xlim([-0.5 10.5]);
ylim([0 2000]); 
title('R^2 = 0.01');
set(gca,'YTick',[]);

set(findall(gcf,'type','axes'),'fontsize',12);
set(findall(gcf,'type','text'),'fontSize',12);
set(findall(gcf,'type','legend'),'fontSize',12);
axesHandles = findobj(get(hist_sel_large,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print('hist_sel_large','-dpdf','-bestfit');
close(hist_sel_large);
