function K = compute_kernel(X, kernel)


%kernel = 'poly';
gamma = 1;
degree = 2;

%defParam

nm = @(X,p)repmat(sum(X.^2,2),1,p);
linKer = @(X1,X2)X1*X2';
rbfKer = @(X1,X2)exp(-(nm(X1,size(X2,1))+nm(X2,size(X1,1))'-2*X1*X2')/2/gamma^2);
lapKer = @(X1,X2)exp(-pdist2(X1,X2)/gamma);
polyKer = @(X1,X2)(1+gamma*X1*X2').^degree;

if strcmpi(kernel,'lin'), kerFun = linKer;
elseif strcmpi(kernel,'poly'), kerFun = polyKer;
elseif strcmpi(kernel,'rbf'), kerFun = rbfKer;
elseif strcmpi(kernel,'lap'), kerFun = lapKer; % Laplacian
else error('unknown ker');
end

K = kerFun(X,X);