function data = rfspectra(images,rf,varargin)
%% RFSPECTRA takes an image series and plots


function img = rfboxplot(dout,rf)
img = zeros([341,1]);
slicewidth = 100;
for i=1:length(dout)
    slice = crop_custom(dout(i).img,[200 80],[340,slicewidth]);
    slice_avg = mean(slice,2);
    img = [img,slice_avg];
end
img = img(:,2:end);
imagesc(img);
set(gca,'XTickLabel',num2str([rf(1:2:end)]'));
end

