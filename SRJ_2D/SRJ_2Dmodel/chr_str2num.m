function chr_num = chr_str2num(chr_str0)

for c1=1:length(chr_str0)

    chr_str=char(chr_str0(c1));
    chri='';
    if length(chr_str)<3,
        chri=chr_str;
    elseif strcmp(chr_str(1:3),'chr')
        chri=chr_str(4:end);
    else
        chr_num(c1)=-1;
    end
    if strcmp(chri,'X')
        chr_num(c1)=23;
    elseif strcmp(chri,'Y')
        chr_num(c1)=24;   
    elseif ~isnan(str2double(chri)),
        chr_num(c1)=str2double(chri);
    else
        chr_num(c1)=-1;
    end

end