function images = listfits(directory)
%% LISTFITS returns a cell array of full paths to fits files in a given directory specified in the
% directory string

x=dir(directory);
names= {x.name};
[fitslocs,~]=cellfun(@size,strfind(names,'fits'));
fitsfiles=x(logical(fitslocs));
images = {fitsfiles.name};
images = cellfun(@(x) [directory,'/',x],images,'UniformOutput',false); % append full path 

end