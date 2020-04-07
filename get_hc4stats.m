function [t,pval,B,V] = get_hc4stats(r,zr,zX,zY)

% sub-routine for compute the t,pval and CI using the HC4 estimator of the
% variance 

[n,p]=size(zX);
if size(r) == [1 1]
    r  = repmat(r,1,p);
    zr = repmat(zr,1,p);
end

for column = p:-1:1
    D             = [ones(n,1) zX(:,column)];
    Betas         = pinv(D)*zY(:,column);
    B(column)     = Betas(2);
    residuals     = zY(:,column) - D*Betas;
    for row=n:-1:1
        h(row,:)  = D(row,:)*inv(D'*D)*D(row,:)';
    end
    d             = min(4,h/mean(h)); %n*h / sum(h)
    S             = inv(D'*D)*D'*diag((residuals.^2)./((1-h).^d))*D*inv(D'*D);
    S             = diag(S); 
    V(column)     = S(2:end); % estimates of squared standard error of r
    t(column)     = r(column)/sqrt(V(column));
    pval(column)  = 2*tcdf(-abs(t(column)),n-2);
end