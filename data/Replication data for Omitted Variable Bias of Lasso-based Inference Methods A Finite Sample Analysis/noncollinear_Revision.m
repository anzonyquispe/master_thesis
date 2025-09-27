function [xN,keep] = noncollinear_Revision(X)

[~,Rx,ex] = qr(X,'vector');
keep = ex(abs(diag(Rx)) > 1e-8);
xN = X(:,keep);
% Make sure not perfectly collinear with a constant
if size(xN,2)+1 > rank([xN ones(size(xN,1),1)])
    ctemp = xN\ones(size(xN,1),1);
    xN(:,ctemp == max(ctemp)) = [];
    keep(ctemp == max(ctemp)) = [];
end
