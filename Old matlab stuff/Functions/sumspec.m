function imgsum = sumspec(images)
%% RFLOAD loads the raw images as OD arrays
    % Initialize data struct
    imgsum = fitsread(images{1});
    % Load the images from the filenames
    fprintf('\n');
%     s2img = loadfitsimage('/Users/biswaroopmukherjee/Documents/Physics/Research/Zwierlein/box data/11-25-2015_19_52_49_top.fits');
    for i =2:length(images)
        fprintf('.');
        imgsum=imgsum +fitsread(images{i});
    end
    fprintf('\n');
end
