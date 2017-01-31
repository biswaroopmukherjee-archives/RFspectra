function sumimg = imagesum(data)
sumimg = data(1).img;
for i=2:length(data)
    sumimg = sumimg + data(i).img;
end