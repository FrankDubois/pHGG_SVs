function TbyGene=TophitsByGenes(hitstable,hitstable_lookup,pos_pad,bins,refgene,refgene_tble,UTumor,CosmicCencus,uFusionTable,bins_annot)

bin_size_th=2e6;
max_num_gene=2; % max number of genes to list per tile

    
TbyGene = struct('bins',[],'gene_i',[],'pos_i',[],'chr_i',[],'pos_i_s',[],'pos_i_e',[],'gene_j',[],'pos_j',[],'chr_j',[],'pos_j_s',[],'pos_j_e',[],'p_val',[],'num_events',[],'annot',[],'avg_dist',[],'std_i',[],'std_j',[],'tumor_1',[],'tumor_2',[],'tumor_3',[],'tumor_1n',[],'tumor_2n',[],'tumor_3n',[],'pp',[],'pm',[],'mp',[],'mm',[]);
isemptyUcode = ~isempty(UTumor);
if isemptyUcode
    UcodeArray=table2array(UTumor);
    [tier1,tier2]=load_tumor_subtype(UTumor);
    [u_tier2,ia,ic_tier2]=unique(tier2);
end

firstline=1;
lastline=0;
Tc=1;

counted_bins = ones(size(hitstable_lookup,1),1);

while sum(counted_bins)>0,
    
    current_bin=[];
    current_bin_idx=find(counted_bins,1);
    counted_bins(current_bin_idx)=0;
    current_bin=hitstable_lookup(current_bin_idx,1:2);
    clines=hitstable_lookup(current_bin_idx,3):hitstable_lookup(current_bin_idx,4);
    TbyGene(Tc).bins(1)=current_bin_idx;
    
    % find adjacent bins with significant events
    if current_bin(1)-1 >0 && current_bin(1)+1 < length(bins)
        contig_bins=counted_bins&hitstable_lookup(:,1)==current_bin(1)+1&hitstable_lookup(:,2)==current_bin(2)&bins(current_bin(1,1))==bins(current_bin(1)+1,1);
        contig_bins=contig_bins|(counted_bins&hitstable_lookup(:,1)==current_bin(1)-1&hitstable_lookup(:,2)==current_bin(2)&bins(current_bin(1,1))==bins(current_bin(1)+1,1));
        contig_bins=contig_bins|(counted_bins&hitstable_lookup(:,1)==current_bin(1)&hitstable_lookup(:,2)==current_bin(2)+1&bins(current_bin(1,1))==bins(current_bin(1)+1,1));
        contig_bins=contig_bins|(counted_bins&hitstable_lookup(:,1)==current_bin(1)&hitstable_lookup(:,2)==current_bin(2)-1&bins(current_bin(1,1))==bins(current_bin(1)+1,1));
        while sum(contig_bins)>0
            sc=sum(contig_bins);
            nnz_bins=find(contig_bins);
            TbyGene(Tc).bins(length(TbyGene(Tc).bins)+1:length(TbyGene(Tc).bins)+length(nnz_bins))=nnz_bins;
            for c1=1:sum(contig_bins),
                clines=[clines hitstable_lookup(nnz_bins(c1),3):hitstable_lookup(nnz_bins(c1),4)];
            end
            counted_bins(contig_bins)=0;
            contig_bins=0;
            for c1=1:sc,           
                current_bin_n=hitstable_lookup(nnz_bins(c1),1:2);
                contig_bins=contig_bins|(counted_bins&hitstable_lookup(:,1)==current_bin_n(1)+1&hitstable_lookup(:,2)==current_bin_n(2)&bins(current_bin_n(1,1))==bins(current_bin_n(1)+1,1));
                contig_bins=contig_bins|(counted_bins&hitstable_lookup(:,1)==current_bin_n(1)-1&hitstable_lookup(:,2)==current_bin_n(2)&bins(current_bin_n(1,1))==bins(current_bin_n(1)+1,1));
                contig_bins=contig_bins|(counted_bins&hitstable_lookup(:,1)==current_bin_n(1)&hitstable_lookup(:,2)==current_bin_n(2)+1&bins(current_bin_n(1,1))==bins(current_bin_n(1)+1,1));
                contig_bins=contig_bins|(counted_bins&hitstable_lookup(:,1)==current_bin_n(1)&hitstable_lookup(:,2)==current_bin_n(2)-1&bins(current_bin_n(1,1))==bins(current_bin_n(1)+1,1));
            end
        end
    end
    %%%%%%%%%%%%%
    
    gene_i={}; gene_j={};
    gene_i_pos=[]; gene_j_pos=[];
    range(Tc,1:3)=[hitstable(clines(1),8) min(hitstable(clines,9)) max(hitstable(clines,9))];          
    [gene_i,gene_i_pos]=get_gene_list(range(Tc,1:3),refgene_tble,refgene,pos_pad,pos_pad);
    [gene_locus, gene_i_pos0]=gene_locus_search(range(Tc,1:3),pos_pad,pos_pad);
    if ~isempty(gene_locus),
        gene_i=gene_locus;
        gene_i_pos=gene_i_pos0(2:3);
    end
        
    range(Tc,4:6)=[hitstable(clines(1),10) min(hitstable(clines,11)) max(hitstable(clines,11))];  
    [gene_j,gene_j_pos]=get_gene_list(range(Tc,4:6),refgene_tble,refgene,pos_pad,pos_pad);
    [gene_locus, gene_j_pos0]=gene_locus_search(range(Tc,4:6),pos_pad,pos_pad);
    if ~isempty(gene_locus),
        gene_j=gene_locus;
        gene_j_pos=gene_j_pos0(2:3);
    end
    
    % find other bins included in the above gene list
    if ~isempty(gene_i_pos) & ~isempty(gene_j_pos)
        i_bins0=find(bins(:,1)==range(Tc,1)&bins(:,3)>gene_i_pos(1)&bins(:,2)<gene_i_pos(2));
        j_bins0=find(bins(:,1)==range(Tc,4)&bins(:,3)>gene_j_pos(1)&bins(:,2)<gene_j_pos(2));
        if min(i_bins0)<min(j_bins0),
            i_bins=i_bins0;
            j_bins=j_bins0;
        else
            i_bins=j_bins0;
            j_bins=i_bins0;
        end
        if ~isempty(i_bins) && ~isempty(j_bins)
            add_tiles = counted_bins&hitstable_lookup(:,1)>=i_bins(1)&hitstable_lookup(:,1)<=i_bins(end)&hitstable_lookup(:,2)>=j_bins(1)&hitstable_lookup(:,2)<=j_bins(end);
            if sum(add_tiles)>0,
                nnz_bins=find(add_tiles);
                TbyGene(Tc).bins(length(TbyGene(Tc).bins)+1:length(TbyGene(Tc).bins)+length(nnz_bins))=nnz_bins;
                for c1=1:sum(add_tiles),
                    clines=[clines hitstable_lookup(nnz_bins(c1),3):hitstable_lookup(nnz_bins(c1),4)];
                end
                counted_bins(add_tiles)=0;
                range(Tc,1:3)=[hitstable(clines(1),8) min(hitstable(clines,9)) max(hitstable(clines,9))];
                range(Tc,4:6)=[hitstable(clines(1),10) min(hitstable(clines,11)) max(hitstable(clines,11))]; 
            end
        end
    end

    [gene_i, gene_j] = known_cancer_genes_annot(gene_i, gene_j, uFusionTable, CosmicCencus);
    
    TbyGene(Tc).pos_i = sprintf('%2d:%9d-%9d',range(Tc,1),range(Tc,2),range(Tc,3));
    TbyGene(Tc).pos_j = sprintf('%2d:%9d-%9d',range(Tc,4),range(Tc,5),range(Tc,6));
    TbyGene(Tc).chr_i=range(Tc,1);
    TbyGene(Tc).pos_i_s=range(Tc,2);
    TbyGene(Tc).pos_i_e=range(Tc,3);
    TbyGene(Tc).chr_j=range(Tc,4);
    TbyGene(Tc).pos_j_s=range(Tc,5);
    TbyGene(Tc).pos_j_e=range(Tc,6);
    if ~isempty(gene_i),
        if length(gene_i)<=max_num_gene,
            TbyGene(Tc).gene_i = gene_i;
        else
            TbyGene(Tc).gene_i = gene_i(1:max_num_gene);
        end
    end
    if ~isempty(gene_j),
        if length(gene_j)<=max_num_gene,
            TbyGene(Tc).gene_j = gene_j;
        else
            TbyGene(Tc).gene_j = gene_j(1:max_num_gene);
        end
    end
    
    if ~isempty(bins_annot),
        if sum(bins_annot(current_bin,2))>0
            TbyGene(Tc).annot='L1';
        elseif sum(bins_annot(current_bin,3))>0
            TbyGene(Tc).annot='FG';
        elseif sum(bins_annot(current_bin,6))>0
            TbyGene(Tc).annot='BL';
        end
    end
       

    
%     add known cancer genes annotation
%     known_genes_i=zeros(length(TbyGene(Tc).gene_i),2);
%     known_genes_j=zeros(length(TbyGene(Tc).gene_j),2);
%     for c1=1:length(TbyGene(Tc).gene_i),
%         for c2=1:length(TbyGene(Tc).gene_j),
%             if (sum(strcmp(TbyGene(Tc).gene_i(c1),uFusionTable(:,1)))>0&sum(strcmp(TbyGene(Tc).gene_j(c2),uFusionTable(:,2)))>0)||(sum(strcmp(TbyGene(Tc).gene_j(c2),uFusionTable(:,1)))>0&sum(strcmp(TbyGene(Tc).gene_i(c1),uFusionTable(:,2)))>0)
%                 known_genes_i(c1,1)=1;
%                 known_genes_j(c2,1)=1;
%             elseif sum(strcmp(TbyGene(Tc).gene_i(c1),CosmicCencus(:,1)))>0 || sum(strcmp(TbyGene(Tc).gene_j(c2),CosmicCencus(:,1)))>0
%                 if sum(strcmp(TbyGene(Tc).gene_i(c1),CosmicCencus(:,1)))>0
%                     known_genes_i(c1,2)=1;
%                 end
%                 if sum(strcmp(TbyGene(Tc).gene_j(c2),CosmicCencus(:,1)))>0 
%                     known_genes_i(c1,2)=1;
%                 end
%             end
%         end
%     end
%     for c1=1:length(TbyGene(Tc).gene_i),
%         if known_genes_i(c1,1)==1,
%             TbyGene(Tc).gene_i(c1)=strcat([TbyGene(Tc).gene_i(c1) '**']);
%         elseif known_genes_i(c1,2)==1,
%             TbyGene(Tc).gene_i(c1)=strcat([TbyGene(Tc).gene_i(c1) '*']);
%         end
%     end
%     for c1=1:length(TbyGene(Tc).gene_j),
%         if known_genes_j(c1,1)==1,
%             TbyGene(Tc).gene_j(c1)=strcat([TbyGene(Tc).gene_j(c1) '**']);
%         elseif known_genes_j(c1,2)==1,
%             TbyGene(Tc).gene_j(c1)=strcat([TbyGene(Tc).gene_j(c1) '*']);
%         end
%     end
%     %%%%%%%%%%%%%%%

    
    TbyGene(Tc).num_events = length(unique(hitstable(clines,9)));
    TbyGene(Tc).p_val = sprintf('%1.1d',min(hitstable(clines,4)));
    TbyGene(Tc).q_val = sprintf('%1.1d',min(hitstable(clines,5)));
    if range(Tc,1)==range(Tc,4),
        TbyGene(Tc).avg_dist = mean(abs((hitstable(clines,9)-hitstable(clines,11))));
    else
        TbyGene(Tc).avg_dist = -1;
    end
    TbyGene(Tc).std_i = std(hitstable(clines,9));
    TbyGene(Tc).std_j = std(hitstable(clines,11));
    if isemptyUcode
        %[tumor_hist_val,tumor_hist_loc]=sort(histc(hitstable(clines,16),1:length(UcodeArray)),'descend');
        [tumor_hist_val,tumor_hist_loc]=sort(histc(ic_tier2(hitstable(clines,16)),1:length(u_tier2)),'descend');
        TbyGene(Tc).tumor_1 = strcat(u_tier2(tumor_hist_loc(1)),' (',num2str(tumor_hist_val(1)),')');
        TbyGene(Tc).tumor_2 = strcat(u_tier2(tumor_hist_loc(2)),' (',num2str(tumor_hist_val(2)),')');
        TbyGene(Tc).tumor_3 = strcat(u_tier2(tumor_hist_loc(3)),' (',num2str(tumor_hist_val(3)),')');
        TbyGene(Tc).tumor_1n = tumor_hist_val(1);
        TbyGene(Tc).tumor_2n = tumor_hist_val(2);
        TbyGene(Tc).tumor_3n = tumor_hist_val(3);
    end
    TbyGene(Tc).pp = sum(((hitstable(clines,14)-1)*2+hitstable(clines,15))==1);
    TbyGene(Tc).mp = sum(((hitstable(clines,14)-1)*2+hitstable(clines,15))==2);
    TbyGene(Tc).pm = sum(((hitstable(clines,14)-1)*2+hitstable(clines,15))==3);
    TbyGene(Tc).mm = sum(((hitstable(clines,14)-1)*2+hitstable(clines,15))==4);
    Tc=Tc+1;

end

sum(counted_bins)
sum([TbyGene.num_events])
%TbySz_tbls=struct2table(TbyGene);
%sort_TbySz_tbls=sortrows(TbySz_tbls,'num_events','descend');

%struct2table(T)
%writetable(struct2table(T),'hitsgenetable.csv','Delimiter',',','QuoteStrings',true);
    
    