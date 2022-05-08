function [gene_symb, gene_pos]=gene_locus_search(loc,pad_low,pad_high)

gene_symb={};
gene_pos=[];
IGH_pos=[14 105566277 106879844];
IGL_pos=[22 22026076 22922913];
IGK_pos=[2 88857361 90235368];
    
if loc(1)==IGH_pos(1) & loc(3)+pad_high>IGH_pos(2) & loc(2)-pad_low<IGH_pos(3)
    gene_symb={'IGH'};
    gene_pos=IGH_pos;
end

if loc(1)==IGL_pos(1) & loc(3)+pad_high>IGL_pos(2) & loc(2)-pad_low<IGL_pos(3)
    gene_symb={'IGL'};
    gene_pos=IGL_pos;
end

if loc(1)==IGK_pos(1) & loc(3)+pad_high>IGK_pos(2) & loc(2)-pad_low<IGK_pos(3)
    gene_symb={'IGK'};
    gene_pos=IGK_pos;
end