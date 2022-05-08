function [hitstable,hitstable_lookup] = HitsTableCV(mfull,pa, pval_tophits, bins_event_tble, qFDR, events, refgene_tble)

gene_pad_neg = 1e4;
gene_pad_pos = 1e4;

disp('generating numeric hits table...');

mat_size=size(pa(:,:,1));
bin_event_matrix=[bins_event_tble events]; % not a matrix but a masterlist 
if iscell(mfull)
    mfull0=mfull{1}+mfull{2}+mfull{3}+mfull{4};
else
    mfull0=mfull;
end
pa0=sum(pa,3);

ct=1;
ct_last=1;
for c1=1:length(pval_tophits),2
    numhits=mfull0(pval_tophits(c1,2));
    event_loc=find(bin_event_matrix(:,3)==pval_tophits(c1,2));
    if numhits~=length(event_loc),
        disp(['wrong number of hits at: ' num2str(pval_tophits(c1,2)) '- -' num2str(numhits) '- -' num2str(length(event_loc))]);
    end
    for c2=1:length(event_loc),
    [bini,binj]=ind2sub(mat_size,pval_tophits(c1,2));    
    hitstable(ct,1:2)=[bini,binj]; %bini/j
    hitstable(ct,3)=numhits; %number of events in bin EQVUIALENT to the number of samples because we only allow 1 sample per bin (also checked this using data.table and the sample # in R) 
    hitstable(ct,4)=pval_tophits(c1,1); %p-value
    hitstable(ct,5)=qFDR(c1); %q-value
    hitstable(ct,6)=pa0(bini,binj); %bin probability
    hitstable(ct,7)=c2; %event #
    hitstable(ct,8)=bin_event_matrix(event_loc(c2),4); %chr i
    hitstable(ct,9)=bin_event_matrix(event_loc(c2),5); %pos i
    hitstable(ct,10)=bin_event_matrix(event_loc(c2),7); %chr j
    hitstable(ct,11)=bin_event_matrix(event_loc(c2),8); %pos j    
    hitstable(ct,12)=getlocusid(hitstable(ct,8),hitstable(ct,9),refgene_tble,gene_pad_neg,gene_pad_pos); %locus_id i
    hitstable(ct,13)=getlocusid(hitstable(ct,10),hitstable(ct,11),refgene_tble,gene_pad_neg,gene_pad_pos); %locus_id j
    hitstable(ct,14)=bin_event_matrix(event_loc(c2),6); %strand i
    hitstable(ct,15)=bin_event_matrix(event_loc(c2),9); %strand j
    hitstable(ct,16)=bin_event_matrix(event_loc(c2),10); %tumor type
    hitstable(ct,17)=bin_event_matrix(event_loc(c2),13); %patient #
    hitstable(ct,18)=bin_event_matrix(event_loc(c2),11); %event #
    hitstable(ct,19)=bin_event_matrix(event_loc(c2),12); %sample #
    ct=ct+1;
    end
    hitstable_lookup(c1,1:2)=[bini,binj];
    hitstable_lookup(c1,3:4)=[ct_last,ct-1];
    ct_last=ct;
end
%save(strcat('hitstable_',num2str(date)),'hitstable');

return
