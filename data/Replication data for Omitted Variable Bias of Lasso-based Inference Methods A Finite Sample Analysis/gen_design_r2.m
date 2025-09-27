function [y, d, X] = gen_design_r2(n,p,k,RsqY,RsqD)

tbeta = zeros(p,1);
tbeta(1:k) = 1;

tgamma = zeros(p,1);
tgamma(1:k) = 1;

cD = sqrt(RsqD/(k-RsqD*k));
cY = sqrt(RsqY/(k-RsqY*k));

gamma = cD*tgamma;
beta = cY*tbeta;

X = randn(n,p);
d = X*gamma + randn(n,1);
y = X*beta + randn(n,1);

end