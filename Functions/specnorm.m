function specout = specnorm(spec)
s=size(spec);
for i=1:s(2)
    specout(:,i) = spec(:,i)/sum(spec(:,i));
end
end
