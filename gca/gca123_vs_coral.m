
% Main routine for comparing Geodesic Covariance Alignment (GCA) with
% correlational alignment (CORAL) 
% Sample invocation of this file: gca123_vs_coral('dslr', 'amazon', 20, 0.9, 0.2, 0.1)

function [accy_coral_mda] = gca123_vs_coral(src, tgt, round, t, w1, w2)

% add path to SVM classifer using LIBSVM 
addpath('../libsvm-3.20/matlab');

% load in Office dataset 
%--------------------I. prepare data--------------------------------------
load(['../data/' src '_SURF_L10.mat']);     % source domain
fts = fts ./ repmat(sum(fts,2),1,size(fts,2));
Source = zscore(fts,1);    clear fts
Source_lbl = labels;           clear labels

load(['../data/' tgt '_SURF_L10.mat']);     % target domain
fts = fts ./ repmat(sum(fts,2),1,size(fts,2));
Target = zscore(fts,1);     clear fts
Target_lbl = labels;            clear labels

Source = double(Source);
Target = double(Target);

fprintf('\nsource (%s) --> target (%s):\n', src, tgt);
fprintf('round     accuracy\n');
%--------------------II. run experiments----------------------------------

nPerClass = 20;

% for fast nearest neighbor search using Open TS Tools package 
% this assumes the MEX files have been compiled for your architecture 
cOpts.UseFastNNSearch = 1; 

% parallel iterations over trials -- needs Parallel Distributed Computing
% Toolbox 
parfor iter = 1 : round
    fprintf('%4d', iter);
    
    inds = split(Source_lbl, nPerClass);
    
    Xr = Source(inds,:);
    Yr = Source_lbl(inds);
    

    Xtt = Target;
    Ytt = Target_lbl;
    
    n = size(Xr, 1); 
    m = size(Xtt,1); 
    
    % compute graph Laplacian kernels for source and target domains
    [vGraph_s, vEigenVecsSub_s, vEigenVals_s, vOpts, vLaplacian_s, vNNInfo_s,vDiffusion_s,vQ_s] = FastLaplacianDetEigs( Xr, cOpts );
    [vGraph_t, vEigenVecsSub_t, vEigenVals_t, vOpts, vLaplacian_t, vNNInfo_t,vDiffusion_t,vQ_t] = FastLaplacianDetEigs( Xtt, cOpts );
    
    % compute source and target graph kernels using manifold diffusion for
    % now
    kernel_s = zeros(size(vEigenVecsSub_s,2)); 
    
    for i=1:size(vEigenVecsSub_s,1) 
        kernel_s = kernel_s + exp(-1/2*vEigenVals_s(i))*vEigenVecsSub_s(i,:)'*vEigenVecsSub_s(i,:); 
    end
    
    kernel_t = zeros(size(vEigenVecsSub_t,2)); 
    
    for i=1:size(vEigenVecsSub_t,1) 
        kernel_t = kernel_t + exp(-1/2*vEigenVals_t(i))*vEigenVecsSub_t(i,:)'*vEigenVecsSub_t(i,:); 
    end
    
    invK = zeros(n+m); 
    
    invK(1:n, 1:n) = inv(kernel_s + eye(n)); % to ensure SPD 
    invK(n+1:m+n, n+1:n+m) = inv(kernel_t + eye(m)); 
    
    
    cov_source = cov(Xr) + eye(size(Xr,2)); 
    cov_target = cov(Xtt) + eye(size(Xtt,2)); 
    
    % form MMD L matrix to minimize distributional differences between
    % source and target
    
    Xall = [Xr; Xtt]; 
    
    Lmat = ones(m + n)*(-1/n*m); 
    
    Lmat(1:n,1:n) = ones(n,n)*(1/n^2); 
    Lmat(n+1:m+n,n+1:m+n) = ones(m,m)*(1/m^2); 
    
    
    A_gca_inv_source = sharp(inv(cov_source + Xall'*(w1*Lmat + w2*invK)*Xall),cov_target, t); 
    
    % Regular CORAL method
    A_coral = cov_source^(-1/2)*cov_target^(1/2);
    
    Sim_coral = Xr * A_coral * Xtt';
    accy_coral(iter) = SVM_Accuracy(Xr, A_coral, Ytt, Sim_coral, Yr);
    
    
    Sim_gca_inv_source = Xr * A_gca_inv_source * Xtt';
    accy_coral_mda(iter) = SVM_Accuracy(Xr, A_gca_inv_source, Ytt, Sim_gca_inv_source, Yr);
    
end

figure; plot(accy_coral_mda, 'r-', 'LineWidth', 2); hold on; plot(accy_coral, 'b-.', 'LineWidth', 2); legend('GCA', 'CORAL'); grid on; 
str = ['Transfer from ', src, ' to target ', tgt]; title(str, 'Fontsize', 16); 

end

function res = SVM_Accuracy (trainset, M,testlabelsref,Sim,trainlabels)
Sim_Trn = trainset * M *  trainset';
index = [1:1:size(Sim,1)]';
Sim = [[1:1:size(Sim,2)]' Sim'];
Sim_Trn = [index Sim_Trn ];

C = [0.001 0.01 0.1 1.0 10 100 1000 10000];
parfor i = 1 :size(C,2)
    model(i) = svmtrain(trainlabels, Sim_Trn, sprintf('-t 4 -c %d -v 2 -q',C(i)));
end
[val indx]=max(model);
CVal = C(indx);
model = svmtrain(trainlabels, Sim_Trn, sprintf('-t 4 -c %d -q',CVal));
[predicted_label, accuracy, decision_values] = svmpredict(testlabelsref, Sim, model);
res = accuracy(1,1);
end


function acc = LinAccuracy(trainset,testset,trainlbl,testlbl)
model = trainSVM_Model(trainset,trainlbl);
[predicted_label, accuracy, decision_values] = svmpredict(testlbl, testset, model);
acc = accuracy(1,1);
end

function svmmodel = trainSVM_Model(trainset,trainlbl)
C = [0.001 0.01 0.1 1.0 10 100 ];
parfor i = 1 :size(C,2)
    model(i) = svmtrain(double(trainlbl), sparse(double((trainset))),sprintf('-c %d -q -v 2',C(i) ));
end
[val indx]=max(model);
CVal = C(indx);
svmmodel = svmtrain(double(trainlbl), sparse(double((trainset))),sprintf('-c %d -q',CVal));
end

function [idx1 idx2] = split(Y,nPerClass, ratio)
% [idx1 idx2] = split(X,Y,nPerClass)
idx1 = [];  idx2 = [];
for C = 1 : max(Y)
    idx = find(Y == C);
    rn = randperm(length(idx));
    if exist('ratio')
        nPerClass = floor(length(idx)*ratio);
    end
    idx1 = [idx1; idx( rn(1:min(nPerClass,length(idx))) ) ];
    idx2 = [idx2; idx( rn(min(nPerClass,length(idx))+1:end) ) ];
end
end
