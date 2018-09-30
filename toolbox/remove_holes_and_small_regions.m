function region_mod = remove_holes_and_small_regions(region, size, background)

if nargin <3
background = 'weak';
end

seg = region;
seg = bwareaopen(seg, size);
seg = 1 - seg;
if strcmp( background,'strong')
    L = bwlabel( seg );
    label = unique(L);
    for i = 1:length(label)
        n(i) = sum(sum(L==label(i)));
    end
    [~, b] = max(n);
    n(b)=0;
    [~, c] = max(n);
    valid = logical((L==label(b))+(L==label(c)));
    seg = seg.*valid;
else
    seg = bwareaopen(seg, size);
end
    region_mod = 1-seg;
end