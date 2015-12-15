function [imgout,data] = sumimg(images,varargin)
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
    case 1 % load samples
        xcrop = 1:512;
        ycrop = 1:512;
    case 3 % provide crop coordinates
        crop = varargin{1};
        xcrop = crop(1):crop(2);
        ycrop = crop(3):crop(4);
    otherwise
        msgbox('Check your parameters');
end



%% Load images
data = rfload(images);

%% Slice images
imgout = rfprocess(data,xcrop,ycrop);



%% Plot spectra (optional)
figure(1)
imagesc(imgout);
axis image
% ax1 = gca;
% set(ax1,'XTick',1:2:length(rf))
% set(ax1,'XTickLabel',num2str(81735-1000*cell2mat(rf(1:2:end)')));
% set(ax1,'FontSize',14);
% xlabel('RF frequency (kHz from 81.735 MHz)');
% ylabel('Axial position');

end


function imout = rfprocess(data,xcrop,ycrop)
%% RFPROCESS rotates and slices the raw images
    % Populate spectra
    imout = zeros(length(xcrop),length(ycrop));
    for i=1:length(data)
        image = imrotate(data(i).img,4); % Rotate the image by 4 degrees
        slice = image(xcrop,ycrop);
        imout = imout + slice;
    end
    
end

function data = rfload(images)
%% RFLOAD loads the raw images as OD arrays
    % Initialize data struct
    % Load the images from the filenames
    fprintf('\n');
%     s2img = loadfitsimage('/Users/biswaroopmukherjee/Documents/Physics/Research/Zwierlein/box data/11-25-2015_19_52_49_top.fits');
    for i =1:length(images)
        fprintf('.');
        data(i).name = images{i};
%         disp(images{i});
%         disp(rf{i});
        data(i).img=loadfitsimage(data(i).name);
    end
    fprintf('\n');
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
