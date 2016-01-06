v = VideoWriter('C:\Users\Elder\Dropbox (MIT)\BEC1\Processed Data\2016\2016-01\2016-01-04\polarized\rf.avi');
v.FrameRate = 5;
open(v);
image = imrotate(data(1).img,4);
   imagesc(image(100:360,150:370));
   axis image;
   colormap parula;
   caxis([-.1 .4]);
set(gca,'nextplot','replacechildren'); 
for k = 1:length(data) 
   image = imrotate(data(k).img,4);
   imagesc(image(100:360,170:370));
   axis image;
   colormap parula;
   titlestring = strcat('\omega - \omega_0 = 2\pi *  ', num2str(1000*(rf{k}-81.736),'%2.1f'),'kHz');
   title(titlestring)
   caxis([-.1 .4]);
   frame = getframe(gcf);
   writeVideo(v,frame);
end

close(v);