function clock = clockplot(img,points,rf)
spec = rfplotter(img,points);
spec_sum = repmat(sum(spec),[20,1]);
spec_norm = spec./spec_sum;
rfrep = repmat(rf,[points,1])';
clock = fliplr(sum(spec_norm.*rfrep));