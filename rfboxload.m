function data = rfboxload
x=dir;
names= {x.name};
[fitslocs,~]=cellfun(@size,strfind(names,'fits'));
fitsfiles=x(logical(fitslocs));
rf =  [81.7100
   81.7180
   81.7220
   81.7260
   81.7300
   81.7320
   81.7340
   81.7360
   81.7380
   81.7420
   81.7460
   81.7500
   81.7550
   81.7600
   81.7700
   81.7800
   81.7900
   81.8000
   81.8100
   81.8200 ];
for i =1:length(fitsfiles)
    disp(i)
    data(i).name = fitsfiles(i).name;
    data(i).img=loadimage(data(i).name);
    data(i).rf = rf(i);
end
end

