function ind = RandomIndex(im)
%This function finds an index on the image given

f = find(im > 0);

i = randi(length(f));

ind = f(i);

end