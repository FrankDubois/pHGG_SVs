function locusid=getlocusid(chr,pos,refgene_lookup,pad_high,pad_low)

loc=find(refgene_lookup(:,1) == chr & refgene_lookup(:,2) < pos & refgene_lookup(:,3) > pos);
if length(loc)>=1;
    locusid=refgene_lookup(loc(1),4);
else
    loc_neg=find(refgene_lookup(:,1) == chr & refgene_lookup(:,3) < pos & refgene_lookup(:,3) > pos - pad_low);
    loc_pos=find(refgene_lookup(:,1) == chr & refgene_lookup(:,2) > pos & refgene_lookup(:,2) < pos + pad_high);
    if isempty(loc_neg) && ~isempty(loc_pos)
        locusid = refgene_lookup(loc_pos(1),4);
    elseif ~isempty(loc_neg) && isempty(loc_pos)
        locusid = refgene_lookup(loc_neg(1),4);
    elseif ~isempty(loc_neg) && ~isempty(loc_pos)
        if abs(refgene_lookup(loc_neg(1),3)-pos) > abs(refgene_lookup(loc_pos(1),2)-pos)
            locusid = refgene_lookup(loc_neg(1),4);
        else
            locusid = refgene_lookup(loc_pos(1),4);
        end
    else
    locusid=-1e10;
    end
end
