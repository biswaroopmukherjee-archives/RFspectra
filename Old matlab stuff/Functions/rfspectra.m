function [data,spec,clocks] = rfspectra(images,rf,numslices,varargin)
%% RFSPECTRA takes an image series and plots
% Usage:  rfspectra(images,rf,crop)
%         images: a cell array with full paths to images
%         rf: a cell array with rf frequencies
%         crop: [x1, x2, y1, y2] where (x1,y1) and (x2,y2) are crop coordinates
%
%         spec: the spectrum output
%         clocks: array of mean RF transition frequencies


%% Arguments
switch nargin
    case 0 % load samples
        [images,rf] = samplesload;
        xcrop = 101:383;
        ycrop = 204:389;
        numslices = 10;
    case 3 % default cropping from 2015-11-18
        xcrop = 100:373;
        ycrop = 264:289;
    case 4 % provide crop coordinates
        crop = varargin{1};
        xcrop = crop(1):crop(2);
        ycrop = crop(3):crop(4);
    otherwise
        msgbox('Check your parameters');
end

rf = reshape(rf,[1 length(rf)]);
images = reshape(images,[1 length(images)]);

% Sort the spectra by RF values
[rfsort,ix] = sort(cell2mat(rf));
images = images(ix);
rf = num2cell(rfsort);

%% Load images
data = rfload(images,rf);

%% Slice images
spec = rfprocess(data,xcrop,ycrop);

%% Bin spectra (optional)
   if numslices~=0
    spec = rfbin(spec,numslices);
   end

%% Find the clock shifts
clocks = clockfind(spec,rf);

%% Plot spectra (optional)
figure(1)
subplot(2,3,2)
imagesc(specnorm(spec));
ax1 = gca;
set(ax1,'XTick',1:2:length(rf))
set(ax1,'XTickLabel',num2str(1000*cell2mat(rf(1:2:end)')-81735));
set(ax1,'FontSize',14);
xlabel('\Delta - 81735 kHz');
ylabel('Axial position');
colormap jet
caxis([0 .2])
title('Axial slices')

%% Plot the clock shifts (optional)
subplot(2,3,3)
plot(clocks,'Marker','.','MarkerSize',15,'LineStyle','none')
ylim([81.735,81.745])
ax2 = gca;
set(ax2,'FontSize',14);
xlabel('Axial position');
ylabel('Mean RF transition frequency');
title('Clock shifts in position')

%% Plot the clock shift as a function of "kf" from spectral summing
specsum = sum(spec,2);
subplot(2,3,1)
plot(specsum, clocks,'Marker','.','MarkerSize',15,'LineStyle','none')
ylim([81.735,81.745])
xlim([0,max(specsum)])
ax3 = gca;
set(ax3,'FontSize',14);
xlabel('k_F (a.u.)');
ylabel('Mean RF transition frequency');
title('Clock shifts in k_F')

%% Plot the cropped image
subplot(2,3,4)
[~,ix]=max(sum(spec));
image = imrotate(data(ix).img,4); % Rotate the image by 4 degrees
imagesc(image(xcrop,ycrop));
axis image
title('Sample cropped image')
    
%% Plot individual spectra
subplot(2,3,5)
plot(spec');
title('Individual spectra')
ax4 = gca;
xlim([1 size(spec,2)])
set(ax4,'XTick',1:2:length(rf))
set(ax4,'XTickLabel',num2str(1000*cell2mat(rf(1:2:end)')-81735));
set(ax4,'FontSize',14);
xlabel('\Delta - 81735 kHz');

%% Plot tails
subplot(2,3,6)
tails = spec(:,end-11:end)';
tailrf = 1000* cell2mat(rf(end-11:end))' - 81735;
plot(tails,'.');
title('tails')
ax4 = gca;
set(ax4,'XTick',1:length(tailrf))
set(ax4,'XTickLabel',num2str(tailrf));
set(ax4,'FontSize',14);
xlabel('\Delta - 81735 kHz');
hold all

for i=1:size(tails,1)
    [xData, yData] = prepareCurveData( tailrf, log(tails(:,i)) );
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft );
    if gof.rsquare>0.5
        plot(exp(fitresult(tailrf)))
    end
end
hold off
end

function clocks = clockfind(spec,rf)
%% CLOCKPLOT plots the clock shift from the spectrum data provided in rf
    s = size(spec);
    specout = specnorm(spec);
    rfrep = repmat(cell2mat(rf),[s(1),1]);
    clocks = sum(specout.*rfrep,2);
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
        spec(:,i) = slice - mean(slice(1:30));
    end
end

function data = rfload(images,rf)
%% RFLOAD loads the raw images as OD arrays
    % Initialize data struct
    data(1:length(images)) = struct('name','','img',[],'rf',0);
    % Load the images from the filenames
    fprintf('\n');
%     s2img = loadfitsimage('C:\2015-12\2015-12-07\hybrid spectrum _1\offres state 2\12-07-2015_20_53_53_top.fits');
    for i =1:length(images)
        fprintf('.');
        data(i).name = images{i};
%         disp(images{i});
%         disp(rf{i});
        data(i).img=loadfitsimage(data(i).name);%-s2img;
        data(i).rf = rf{i};
    end
    fprintf('\n');
end

function [images,rf] = samplesload
%% Load the sample images
    images = {'11-19-2015_01_47_24_top.fits','11-19-2015_01_48_17_top.fits','11-19-2015_01_49_10_top.fits','11-19-2015_01_50_03_top.fits','11-19-2015_01_50_56_top.fits','11-19-2015_01_51_49_top.fits','11-19-2015_01_52_42_top.fits','11-19-2015_01_53_35_top.fits','11-19-2015_01_54_28_top.fits','11-19-2015_01_55_21_top.fits','11-19-2015_01_56_15_top.fits','11-19-2015_01_57_08_top.fits','11-19-2015_01_58_01_top.fits','11-19-2015_01_58_54_top.fits','11-19-2015_01_59_47_top.fits','11-19-2015_02_00_40_top.fits','11-19-2015_02_01_33_top.fits','11-19-2015_02_03_02_top.fits','11-19-2015_02_03_55_top.fits','11-19-2015_02_04_48_top.fits'};
    directory = '../Samples/';
    images = cellfun(@(x) [directory,x],images,'UniformOutput',false); % append full path 
    rf = num2cell([81.7100000000000;81.7180000000000;81.7220000000000;81.7260000000000;81.7300000000000;81.7320000000000;81.7340000000000;81.7360000000000;81.7380000000000;81.7420000000000;81.7460000000000;81.7500000000000;81.7550000000000;81.7600000000000;81.7700000000000;81.7800000000000;81.7900000000000;81.8000000000000;81.8100000000000;81.8200000000000]');
end

function img = loadfitsimage(filename)
 data=fitsread(filename);
    absimg=(data(:,:,2)-data(:,:,3))./(data(:,:,1)-data(:,:,3));

%     % Identify "burned pixels" and make sure they will come out zero.
%     burnedpoints = absimg < 0;
%     absimg(burnedpoints) = 1;
% 
%     % Same thing for points which should be accidentally equal to zero
%     % (withatoms == background) or infinity (withoutatoms == background)
%     zeropoints = absimg == 0;
%     absimg(zeropoints) = 1;
% 
%     infpoints = abs(absimg) == Inf;
%     absimg(infpoints) = 1;
% 
%     nanpoints = isnan(absimg);
%     absimg(nanpoints) = 1;

%replace the pixels with a value of negtive number,0 or inf or nan by the
%average of nearset site.
    ny=size(absimg,1);
    nx=size(absimg,2);
    burnedpoints = absimg <= 0;
    infpoints = abs(absimg) == Inf;
    nanpoints = isnan(absimg);
    Change=or(or(burnedpoints,infpoints),nanpoints);
    NChange=not(Change);
    for i=2:(ny-1)
        for j=2:(nx-1)
            if Change(i,j)
                n=0;
                rp=0;
                if NChange(i-1,j)
                    rp=rp+absimg(i-1,j);
                    n=n+1;
                end
                if NChange(i+1,j)
                    rp=rp+absimg(i+1,j);
                    n=n+1;
                end
                if NChange(i,j-1)
                    rp=rp+absimg(i,j-1);
                    n=n+1;
                end
                if NChange(i,j+1)
                    rp=rp+absimg(i,j+1);
                    n=n+1;
                end
                if (n>0)
                    absimg(i,j)=(rp/n);
                    Change(i,j)=0;
                end
            end
        end
    end
    absimg(Change)=1;
    img = log(absimg);
end
