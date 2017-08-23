function graphannot(charvec, locs,amplitude)
% locs = locations of annotations to add to graph
% charvec = characters that represent labels that will be added to graph 
% amplitude = amplitude to mark labels

if islogical(charvec)
    charvec = +charvec;
    chars = string(charvec);
    chars = char(chars);
elseif ischar(charvec)
    chars = charvec;
elseif isnumeric(charvec)
    chars = string(charvec);
    chars = char(chars(:));
elseif iscell(charvec)
    chars = char(charvec);
end


for i = 1:length(chars)
    txt = chars(i);
    text(locs(i),amplitude,txt);
end

end