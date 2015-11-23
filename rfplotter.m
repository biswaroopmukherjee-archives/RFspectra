function spec = rfplotter(img,slices)
s = size(img);
w = floor(s(1) / slices);
for i=1:slices
    low = (1+(i-1)*w);
    high = i*w;
    spec(:,i) = sum(img(low:high,:))';
end
    