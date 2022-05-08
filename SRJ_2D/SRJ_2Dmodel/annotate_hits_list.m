function annotated_table = annotate_hits_list( TbyGene_Table,SVTable,bins,hitstable_lookup,pa_mix)

if isstruct(TbyGene_Table)
    annotated_table=table();
    chits=table(); ct=1;
    for c1=1:length(TbyGene_Table)
        for c2=1:length(TbyGene_Table(c1).bins)
            clear chits
            hitstable_bin=hitstable_lookup(TbyGene_Table(c1).bins(c2),1:2);
            chri=bins(hitstable_bin(1),1);
            posi=bins(hitstable_bin(1),2:3);
            chrj=bins(hitstable_bin(2),1);
            posj=bins(hitstable_bin(2),2:3);
            % merging bins?
            SVTable_lines = SVTable.seqnames==chri&SVTable.start>=posi(1)&SVTable.start<=posi(2) & ...
                        SVTable.altchr==chrj&SVTable.altpos>=posj(1)&SVTable.altpos<=posj(2) | ...
                        SVTable.altchr==chri&SVTable.altpos>=posi(1)&SVTable.altpos<=posi(2) & ...
                        SVTable.seqnames==chrj&SVTable.start>=posj(1)&SVTable.start<=posj(2);
            chits=SVTable(SVTable_lines,:);
            if isempty(TbyGene_Table(c1).gene_i),
                TbyGene_Table(c1).gene_i = {'none'};
            end
            if isempty(TbyGene_Table(c1).gene_j),
                TbyGene_Table(c1).gene_j = {'none'};
            end
            chits.gene_i = repmat(TbyGene_Table(c1).gene_i(1),height(chits),1);
            chits.gene_j = repmat(TbyGene_Table(c1).gene_j(1),height(chits),1);
            chits.hit_num=c1*ones(sum(SVTable_lines),1);
            chits.pval=str2double(TbyGene_Table(c1).p_val)*ones(sum(SVTable_lines),1);
            chits.qval=str2double(TbyGene_Table(c1).q_val)*ones(sum(SVTable_lines),1);            
            chits.tile_num=c2*ones(sum(SVTable_lines),1);
            chits.u_tile_num=ct*ones(sum(SVTable_lines),1);
            ct=ct+1;
            chits.p_mix=pa_mix(min(hitstable_bin),max(hitstable_bin))*ones(sum(SVTable_lines),1);
            annotated_table=[annotated_table;chits];
        end
    end    
else    
    annotated_table=table();
    chits=table();
    for c1=1:length(TbyGene_Table)        
            clear chits
            hitstable_bin=TbyGene_Table(c1,:);

            chri=bins(hitstable_bin(1),1);
            posi=bins(hitstable_bin(1),2:3);
            chrj=bins(hitstable_bin(2),1);
            posj=bins(hitstable_bin(2),2:3);

            SVTable_lines = SVTable.seqnames==chri&SVTable.start>=posi(1)&SVTable.start<=posi(2) & ...
                            SVTable.altchr==chrj&SVTable.altpos>=posj(1)&SVTable.altpos<=posj(2) | ...
                            SVTable.altchr==chri&SVTable.altpos>=posi(1)&SVTable.altpos<=posi(2) & ...
                            SVTable.seqnames==chrj&SVTable.start>=posj(1)&SVTable.start<=posj(2);

            chits=SVTable(SVTable_lines,:);
            chits.hit_num=c1*ones(sum(SVTable_lines),1);
            chits.pval=str2double(TbyGene_Table(c1).p_val)*ones(sum(SVTable_lines),1);
            chits.num_tiles=size(TbyGene_Table(c1).bins,2)*ones(sum(SVTable_lines),1);
            chits.p_mix=pa_mix(min(hitstable_bin),max(hitstable_bin))*ones(sum(SVTable_lines),1);
            annotated_table=[annotated_table;chits];
        
    end
    
end
