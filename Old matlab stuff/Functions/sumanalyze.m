s2profile = sum(s2image,2);
s2profile = s2profile(50:end)-mean(s2profile(end-50:end));
s2profnorm = s2profile/sum(s2profile);


for i=1:size(data,3)
    image = data(:,:,i);
    profile = sum(image,2);
    profile = profile(50:end)-mean(profile(end-50:end));
    profnorm(:,i) = profile/sum(profile);
    diff(:,i) = profnorm(:,i)-s2profnorm;
end
    