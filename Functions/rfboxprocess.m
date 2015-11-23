function dout = rfboxprocess(data)
for i=1:length(data)
     dout(i).img = crop_custom(imrotate(data(i).img,4),[250,290],[400,150]);
end
end


function imgout=crop_custom(img,r,d)

w=floor(d/2);
x=r(1);
y=r(2);
imgout = img((x-w(1)):(x+w(1)),(y-w(2)):(y+w(2)));
end