function [rg_symbol,pos]=get_gene_list(loci,refgene_lookup,refgene,pad_high,pad_low)

list1=refgene_lookup(refgene_lookup(:,1) == loci(1) & refgene_lookup(:,3) >= loci(2) & refgene_lookup(:,2) <= loci(3),4);

list2=refgene_lookup(refgene_lookup(:,1) == loci(1) & refgene_lookup(:,3) >= loci(2)-pad_low & refgene_lookup(:,2) <= loci(2),4);

list3=refgene_lookup(refgene_lookup(:,1) == loci(1) & refgene_lookup(:,3) >= loci(3) & refgene_lookup(:,2) <= loci(3)+pad_high,4);

list=unique([list1;list2;list3]);

list_prim = list(ismember(list ,list1));
list_sec = list(~ismember(list ,list1));

empty_list_pad = 5e4;
empty_list_pad_th = 5e5;
if isempty(list_prim) && isempty(list_sec)
    while isempty(list_sec) && empty_list_pad<empty_list_pad_th
        list2=refgene_lookup(refgene_lookup(:,1) == loci(1) & refgene_lookup(:,3) >= loci(2)-empty_list_pad & refgene_lookup(:,2) <= loci(3)+empty_list_pad,4);
        list_sec = unique(list2);
        empty_list_pad = empty_list_pad + 5e4;
    end
end
    
rg_symbol={};
pos=[];

for c1=1:length(list_prim),
    refgene_loc=find([refgene.rg(:).locus_id]==list_prim(c1));
    if ~isempty(refgene_loc)
        rg_symbol(c1)={refgene.rg(refgene_loc(1)).symb};
        pos(c1,1)=min([refgene.rg(refgene_loc).start]);
        pos(c1,2)=max([refgene.rg(refgene_loc).end]);
    else
        rg_symbol(c1)={};
        pos(c1)=[];
    end
end

for c1=length(list_prim)+1:length(list_prim)+length(list_sec),
    refgene_loc=find([refgene.rg(:).locus_id]==list_sec(c1-length(list_prim)));
    if ~isempty(refgene_loc)
        rg_symbol(c1)={strcat('(',refgene.rg(refgene_loc(1)).symb,')')};
        pos(c1,1)=min([refgene.rg(refgene_loc).start]);
        pos(c1,2)=max([refgene.rg(refgene_loc).end]);
    else
        rg_symbol(c1)={};
        pos(c1)=[];
    end
end

