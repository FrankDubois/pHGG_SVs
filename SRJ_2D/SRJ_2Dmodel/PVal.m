function [qFDR, pa, pval_tophits,mfull] = PVal(mfull, p, qe, sij,approx_flag)
global FDR_THRESHOLD 
qqplot_flag = 0;
%reproducible results
rng(3014)

% general variables
mat_size = size(mfull);

% keep only upper diagonal part of mfull 
mfull = triu(full(mfull));
%divide all diagonal elements by 2 (because we only want to count events in
%i,i bins once (we added mfull + mfull' which effectively counted them
%twice)
mfull(eye(mat_size)~=0) = diag(mfull)/2;
%
nume = sum(mfull(:));

% create probability matrix
if isempty(p)
    pa = bsxfun(@times,qe,sij);
    pa = triu(pa + pa');
    pa(eye(mat_size)~=0) = diag(pa)/2;
else
    pa=p;
    if issymmetric(p)
        pa=triu(p);
        pa(eye(mat_size)~=0) = diag(pa)/2;
    end
end
    
%export pa and mfull matrices to R from here%
%dlmwrite('/xchip/beroukhimlab/Kiran/adjancencies/20190426mfull.txt', mfull, 'delimiter','\t','newline','pc','precision',13);
%dlmwrite('/xchip/beroukhimlab/Kiran/adjancencies/20190426pa.txt', pa, 'delimiter','\t','newline','pc','precision',13);


% divide tiles with positive values from zeros 
pos_k = find(mfull>=1 & pa>0);
zero_k = find(mfull==0 & pa>0);
mfull_pos = full(mfull(pos_k));
p_pos = pa(pos_k);
p_zero = pa(zero_k);

% p-vals for tiles with zero events
disp(['calculating p-val for ' num2str(length(zero_k)) ' tiles with zero events']);
tic
if ~approx_flag
    pval_low=1-(1-p_zero).^nume.*rand(length(zero_k),1);
else
    pval_low=1-exp(-p_zero*nume).*rand(length(zero_k),1);
end
toc

disp(['calculating p-val for the rest ' num2str(length(pos_k)) ' tiles']);
rand_nnz = rand(length(pos_k),1);
tic
if ~approx_flag
    %adding a random number to the p-value, to smooth the p value distribution for fdr  
    %but makes the bottom of the list stochasatic rather than
    %deterministic
   
    pval_high = binopdf(mfull_pos, nume, p_pos).*rand_nnz+(1-binocdf(mfull_pos,nume,p_pos));
else
    pval_high = poisspdf(mfull_pos, nume*p_pos).*rand_nnz+(1-poisscdf(mfull_pos,nume*p_pos));
end
toc
pval_tble=[pval_low zero_k;pval_high pos_k];


pval=sortrows(pval_tble,1);
len_pval=length(pval);

% qq-plot  
if qqplot_flag
   %qqplot(pval(:,1), ;)
end

% calcualte BH-FDR
qFDR=mafdr(pval(:,1),'BHFDR','true');
%Kiran added FDR_THRESHOLD as a global variable
hits_idx=find(qFDR< FDR_THRESHOLD);
tophits=length(hits_idx);
pval_tophits = pval(1:tophits,:);
 
return

