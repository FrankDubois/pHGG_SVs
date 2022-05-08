%set directory path% 
pwd = '/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/matlab/SRJ_2Dmodel/' 
addpath(genpath(pwd));
%% 

%global variable for FDR threshold%
global FDR_THRESHOLD
FDR_THRESHOLD = 0.25;
%% 

load ICGC_2D_SV_model.mat
%% 

sv_file='/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200207_2dmodel/20200207_SV.txt';

SVTable=readtable(sv_file);
EventLengthThreshold=1e3;
%% 

% events list
events0=zeros(height(SVTable),12);

if isa(SVTable.seqnames,'numeric')    
events0(:,1)=SVTable.seqnames;
 else
   events0(:,1)=chr_str2num(SVTable.seqnames)';
   SVTable.seqnames=events0(:,1);
end
%% 

if isa(SVTable.altchr,'numeric')
    events0(:,4)=SVTable.altchr;    
else
    events0(:,4)=chr_str2num(SVTable.altchr)';
SVTable.altchr=events0(:,4);
end
%% 

events0(:,2)=SVTable.start;
[Ustrand1, ia_strand1, ic_strand1]=unique(SVTable.strand);
events0(:,3)=ic_strand1;

events0(:,5)=SVTable.altpos;
[Ustrand2, ia_strand2, ic_strand2]=unique(SVTable.altstrand);
events0(:,6)=ic_strand2;

%[UTumor, ia_code, ic_code]=unique(SVTable.dcc_project_code);
events0(:,7)=1;

[Uevent, ia_event, ic_event]=unique(SVTable.ID);
events0(:,8)=ic_event;

[Usample, ia_sample, ic_sample]=unique(SVTable.sid);
events0(:,9)=ic_sample;
%% 

if sum(strcmp('donor_unique_id',SVTable.Properties.VariableNames))>0
[Upatient, ia_patient, ic_patient]=unique(SVTable.donor_unique_id);
 else
   [Upatient, ia_patient, ic_patient]=unique(SVTable.sid);
end
events0(:,10)=ic_patient;

events0(:,11)=SVTable.HOMLEN;
events0(:,12)=SVTable.INSLEN;

disp(strcat('total events from vcfs: ',num2str(length(events0))));
% filter mask track
[events0,masked_events] = mask_events( events0,mask_track );
disp(strcat('total events after masked regions: ',num2str(length(events0))));
%% 

% events matrix
mfull0=sparse(length(bins),length(bins)); % bins from pan cancer
bins_event_tble0=zeros(length(events0),3);
for c1 = 1:length(events0)
    bini0=find(bins(:,1)==events0(c1,1) & bins(:,2)<=events0(c1,2) & bins(:,3)>=events0(c1,2));  % >sort events into bins 
    binj0=find(bins(:,1)==events0(c1,4) & bins(:,2)<=events0(c1,5) & bins(:,3)>=events0(c1,5)); 
    if ~isempty(bini0) && ~isempty(binj0)
        bins_event_tble0(c1,1) = bini0;
        bins_event_tble0(c1,2) = binj0;
        bini = min(bins_event_tble0(c1,1),bins_event_tble0(c1,2));
        binj = max(bins_event_tble0(c1,1),bins_event_tble0(c1,2));
        bins_event_tble0(c1,3) = sub2ind([length(bins) length(bins)],bini,binj); % note: assign values only to upper tria of the matrix
        mfull0(bini,binj) = mfull0(bini,binj) + 1;
    end
end
%list of pan-cancer bins for each event 
% cols: pan-can bin_i, pan-can bin_j, index in matrix
%% 

% remove unassigned events
unassigned_events=bins_event_tble0(:,1)==0;
events0(unassigned_events,:)=[];
bins_event_tble0(unassigned_events,:)=[];
disp(strcat('total events after removing unassigned events: ',num2str(length(events0))));
%% 

% remove same bin same sample events
T_bin_sample=table();
T_bin_sample.bin=bins_event_tble0(:,3);
T_bin_sample.sample=events0(:,9); %sampleID
[u_t,ia_t,ic_t]=unique(T_bin_sample);
events=events0(ia_t,:);
bins_event_tble=bins_event_tble0(ia_t,:);
disp(strcat('total events after removing same bin same sample events: ',num2str(length(events))));
%% 
events_beforeLengThresh = events;
%EventLengthThreshold = 100;
%% remove below length threshold events
below_length_th=(abs(events(:,5)-events(:,2))<EventLengthThreshold)&(events(:,1)==events(:,4)); % EventLengthThreshold
events(below_length_th,:)=[];
bins_event_tble(below_length_th,:)=[];
disp(strcat('total events after removing below length threshold: ',num2str(length(events))));

%% 
belowLengTh = events_beforeLengThresh(below_length_th,:);

%%
% update mfull
mfull=sparse(length(bins),length(bins));
for c1 = 1:length(events)
    bini0=find(bins(:,1)==events(c1,1) & bins(:,2)<=events(c1,2) & bins(:,3)>=events(c1,2));
    binj0=find(bins(:,1)==events(c1,4) & bins(:,2)<=events(c1,5) & bins(:,3)>=events(c1,5));
    if ~isempty(bini0) && ~isempty(binj0)
        mfull(bini0,binj0) = mfull(bini0,binj0) + 1;
    end
end
%% 

%[qFDR_mix, pa_mix, pval_tophits_mix, mfull_pval_mix] = PVal(mfull, mix_model, bins, events, sij1dx, chsize, CHR, [], [], 1);
[qFDR_mix, pa_mix, FDRpval_tophits, mfull_pval_mix] = PVal(mfull+mfull', mix_model, [], [],1); 
% pa_mix = mix_model = background propability
% mfull_pval_mix = mfull+mfull' = input
%% 

[hitstable_mix,hitstable_mix_lookup] = HitsTableCV(mfull_pval_mix,pa_mix, FDRpval_tophits, bins_event_tble, qFDR_mix, events, refgene_tble);
%%
TbyGene_mix = TophitsByGenes(hitstable_mix,hitstable_mix_lookup,1e4,bins,refgene,refgene_tble,[],CosmicCencus,CuratedFusionGene0,[]);
    
annotated_table = annotate_hits_list( TbyGene_mix,SVTable,bins,hitstable_mix_lookup,pa_mix );
%%
hits_table=table();
hits_table.hit_num = annotated_table.hit_num;
hits_table.gene_i = annotated_table.gene_i;
hits_table.gene_j = annotated_table.gene_j;
hits_table.donor_unique_id = annotated_table.sid;
hits_table.chr_i = annotated_table.seqnames;
hits_table.pos_i = annotated_table.start;
hits_table.strand_i = annotated_table.strand;
hits_table.chr_j = annotated_table.altchr;
hits_table.pos_j = annotated_table.altpos;
hits_table.strand_j = annotated_table.altstrand;
hits_table.pval = annotated_table.pval;
hits_table.uid = annotated_table.uid;
hits_table.EVDNC = annotated_table.EVDNC;
hits_table.QUALITY = annotated_table.QUALITY;
hits_table.svtype = annotated_table.svtype;
writetable(hits_table,'/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200210_2Dmodel/202002102Dhitstbl_q005.txt','delimiter','\t')
writetable(annotated_table,'/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200210_2Dmodel/20200210_2Dannotatedtblq005.txt','delimiter','\t')
 %% 
%  bins_table=table();
%  bins_table.chr = bins(:,1);
%  bins_table.pos1 = bins(:,2);
%  bins_table.pos2 = bins(:,3);
%  bins_table.count = bins(:,4);
%  writetable(bins_table,'/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200210_2Dmodel/20200210_2DpancanBins.txt','delimiter','\t')
%%
% mrfull = full(mfull);
% dlmwrite('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200210_2Dmodel/20200210_2Dtilematrix.txt',mrfull,'delimiter','\t')
dlmwrite('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200210_2Dmodel/20200210_2Dhitstable_mix.txt',hitstable_mix,'delimiter','\t')
%cols(hitstable_mix)= %bini/j, # of samples, p-value, q-value, bin
%probability, chr_i, pos_i, chr_j, pos_j, strand_i, strand_j, tumor_type,
%patient #, event #, 
%%
pa_mix_log = log10(pa_mix +1);
%dlmwrite('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200210_2Dmodel/20200210_pa_mix_log.txt',pa_mix_log,'delimiter','\t')
%dlmwrite('/Users/fdubois/Dropbox (Partners HealthCare)/DIPG/data/20200210_2Dmodel/20200210_pa_mix.txt',pa_mix,'delimiter','\t')
