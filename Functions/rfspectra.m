function [spec,clock] = rfspectra(images,rf,varargin)
%% RFSPECTRA takes an image series and plots
% Usage:  data = rfspectra(images,rf,crop)
%         images: a cell array with full paths to images
%         rf: a cell array with rf frequencies
%         crop: [x1, x2, y1, y2] where (x1,y1) and (x2,y2) are crop coordinates
%
%         specbinned: the spectrum output


%% Arguments
switch nargin
    case 0 % load samples
        [images,rf] = samplesload;
        xcrop = 101:383;
        ycrop = 204:389;
    case 2 % default cropping from 2015-11-18
        xcrop = 111:403;
        ycrop = 294:319;
    case 3 % provide crop coordinates
        crop = varargin{1};
        xcrop = crop(1):crop(2);
        ycrop = crop(3):crop(4);
    otherwise
        msgbox('Check your parameters');
end

%% Load images
data = rfload(images,rf);

%% Slice images
spec = rfprocess(data,xcrop,ycrop);

%% Bin spectra (optional)
% spec = rfbin(spec,30);

%% Find the clock shifts
clock = clockfind(spec,rf);

%% Plot spectra (optional)
figure(1)
imagesc(spec);
ax1 = gca;
set(ax1,'XTick',1:2:length(rf))
set(ax1,'XTickLabel',num2str(cell2mat(rf(2:2:end)),'%.2f'));
set(ax1,'FontSize',14);
xlabel('RF frequency (MHz)');
ylabel('Axial position');

%% Plot the clock shifts (optional)
figure(2);
plot(clock,'Marker','.','MarkerSize',15,'LineStyle','none')
ylim([81.72,81.746])
ax2 = gca;
set(ax2,'FontSize',14);
xlabel('Axial position');
ylabel('Mean RF transition frequency');


end

function clock = clockfind(spec,rf)
%% CLOCKPLOT plots the clock shift from the spectrum data provided in rf
    s = size(spec);
    specout = specnorm(spec);
    rfrep = repmat(cell2mat(rf)',[s(1),1]);
    clock = sum(specout.*rfrep,2);
end

function specout = specnorm(spec)
%% SPECNORM normalizes the spectra
    s=size(spec);
    spec_sum = repmat(sum(spec,2),[1,s(2)]);
    specout = spec./spec_sum;
end

function spec = rfbin(img,slices)
%% RFBIN bins the image along the axial direction
    s = size(img);
    w = floor(s(1) / slices);
    spec = zeros(slices,s(2));
    for i=1:slices
        low = (1+(i-1)*w);
        high = i*w;
        spec(i,:) = mean(img(low:high,:))';
    end
end
    

function spec = rfprocess(data,xcrop,ycrop)
%% RFPROCESS rotates and slices the raw images
    % Initialize spectra
    spec = zeros(length(xcrop),length(data));
    
    % Populate spectra
    for i=1:length(data)
        image = imrotate(data(i).img,4); % Rotate the image by 4 degrees
        slice = mean(image(xcrop,ycrop),2);
        spec(:,i) = slice;
    end
end

function data = rfload(images,rf)
%% RFLOAD loads the raw images as OD arrays
    % Initialize data struct
    data(1:length(images)) = struct('name','','img',[],'rf',0);
    % Load the images from the filenames
    fprintf('\n');
    for i =1:length(images)
        fprintf('.');
        data(i).name = images{i};
        data(i).img=loadimage(data(i).name);
        data(i).rf = rf{i};
    end
    fprintf('\n');
end

function [images,rf] = samplesload
%% Load the sample images
    images = {'11-19-2015_01_47_24_top.fits','11-19-2015_01_48_17_top.fits','11-19-2015_01_49_10_top.fits','11-19-2015_01_50_03_top.fits','11-19-2015_01_50_56_top.fits','11-19-2015_01_51_49_top.fits','11-19-2015_01_52_42_top.fits','11-19-2015_01_53_35_top.fits','11-19-2015_01_54_28_top.fits','11-19-2015_01_55_21_top.fits','11-19-2015_01_56_15_top.fits','11-19-2015_01_57_08_top.fits','11-19-2015_01_58_01_top.fits','11-19-2015_01_58_54_top.fits','11-19-2015_01_59_47_top.fits','11-19-2015_02_00_40_top.fits','11-19-2015_02_01_33_top.fits','11-19-2015_02_03_02_top.fits','11-19-2015_02_03_55_top.fits','11-19-2015_02_04_48_top.fits'};
    directory = '../Samples/';
    images = cellfun(@(x) [directory,x],images,'UniformOutput',false); % append full path 
    rf = num2cell([81.7100000000000;81.7180000000000;81.7220000000000;81.7260000000000;81.7300000000000;81.7320000000000;81.7340000000000;81.7360000000000;81.7380000000000;81.7420000000000;81.7460000000000;81.7500000000000;81.7550000000000;81.7600000000000;81.7700000000000;81.7800000000000;81.7900000000000;81.8000000000000;81.8100000000000;81.8200000000000]);
end

