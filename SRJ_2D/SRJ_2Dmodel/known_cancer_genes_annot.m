function [gene_i, gene_j] = known_cancer_genes_annot(gene_i, gene_j, uFusionTable, CosmicCencus)

len_i=length(gene_i);
len_j=length(gene_j);

if len_i>0 && len_j>0
    known_genes_i=zeros(len_i,2);
    known_genes_j=zeros(len_j,2);
    for c1=1:len_i,
        if strcmp(gene_i{c1}(1),'(')
            gene_i0{c1}=gene_i{c1}(2:end-1);
        else
            gene_i0{c1}=gene_i{c1};
        end
    end
    for c1=1:len_j,
        if strcmp(gene_j{c1}(1),'(')
            gene_j0{c1}=gene_j{c1}(2:end-1);
        else
            gene_j0{c1}=gene_j{c1};
        end
    end
    
    for c1=1:len_i,
        for c2=1:len_j,
            if sum(strcmp(gene_i0(c1),uFusionTable(:,1))&strcmp(gene_j0(c2),uFusionTable(:,2)))>0
                known_genes_i(c1,1)=1;
                known_genes_j(c2,1)=1;
            elseif sum(strcmp(gene_i0(c1),CosmicCencus(:,1)))>0 || sum(strcmp(gene_j0(c2),CosmicCencus(:,1)))>0
                if sum(strcmp(gene_i0(c1),CosmicCencus(:,1)))>0
                    known_genes_i(c1,2)=1;
                end
                if sum(strcmp(gene_j0(c2),CosmicCencus(:,1)))>0 
                    known_genes_j(c1,2)=1;
                end
            end
        end
    end
    for c1=1:len_i,
        if known_genes_i(c1,1)==1,
            gene_i{c1}=strcat([gene_i{c1} '**']);
        elseif known_genes_i(c1,2)==1,
            gene_i{c1}=strcat([gene_i{c1} '*']);
        end
    end
    for c1=1:len_j,
        if known_genes_j(c1,1)==1,
            gene_j{c1}=strcat([gene_j{c1} '**']);
        elseif known_genes_j(c1,2)==1,
            gene_j{c1}=strcat([gene_j{c1} '*']);
        end
    end
end