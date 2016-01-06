s2profile = sum(imrotate(s2image(xcrop,ycrop),4),2);
s2profile = s2profile-mean(s2profile(end-50:end));
s2profnorm = s2profile/sum(s2profile);
image = data(1).img;

for i=2:size(data,2)
    image = image + data(i).img;
end
    
    profile = sum(imrotate(image(xcrop,ycrop),4),2);
    profile = profile-mean(profile(end-50:end));
    profnorm = profile/sum(profile);
    diff = profnorm-s2profnorm;