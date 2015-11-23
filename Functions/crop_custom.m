function imgout=crop_custom(img,r,d)

w=floor(d/2);
x=r(1);
y=r(2);
imgout = img((x-w(1)):(x+w(1)),(y-w(2)):(y+w(2)));
end