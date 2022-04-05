function ind = get_index(description,pattern)
M=length(description);
N=length(pattern);
ind=[];
for i=1:M
    for j=1:N
        if strcmp(description{i},pattern{j})
            ind(end+1)=i;
        end
    end
end
end
