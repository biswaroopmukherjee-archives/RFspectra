function data = RFload(varargin)

%% Define inputs
    switch nargin
        case 0
            exp_name = 'RF_1';
        case 1
            exp_name = varargin{1};
            if ~ischar(exp_name)
                msgbox('Please enter a string for the experiment name');
            end
    end
               
 %% Define paths           
    addpath('Snippet_Readout');

    raw_data_path = 'C:\Users\BEC1\Dropbox (MIT)\BEC1\Image Data and Cicero Files\Data - Raw Images';
    processed_data_path = 'C:\Users\BEC1\Dropbox (MIT)\BEC1\Processed Data';
    
    [sourcedir,tempfolder,destination] = foldermanagement(raw_data_path,processed_data_path,exp_name);

 %% Find files and rf values
 
    [s2file,s2offresfile,s3files,rf] = findfiles(sourcedir)
 
    
 %% Calculate RF spectra
    
    [data,spec,clocks] = rfspectra(s3files,rf);
    
    
    
    
end


function [sourcedir,tempfolder,destination] = foldermanagement(raw_data_path,processed_data_path,exp_name)
% prepare the folders to get raw data from, copy to while working, and
% finally save to after finishing.
    c=clock;
    year = num2str(c(1));
    month = strcat(year,'-',num2str(c(2)));
    day = strcat(month,'-',num2str(c(3),'%02d'));

    sourcedir = strcat(raw_data_path,'\',year,'\',month,'\',day);
    tempfolder = strcat('C:\',year,'\',month,'\',day);
    destination=strcat(processed_data_path,'\',year,'\',month,'\',day,'\',exp_name); % Put the image files in here

    if not(isdir(sourcedir))
        sourcedir = strcat(raw_data_path,'\',year,'\',month);
        if not(isdir(sourcedir))
            sourcedir = strcat(raw_data_path,'\',year);
            if not(isdir(sourcedir))
                sourcedir = strcat(raw_data_path);
            end
        end
    end
    if not(isdir(tempfolder))
        mkdir(tempfolder);
    end
    if not(isdir(destination))
        mkdir(destination);
    end
end

function [s2file,s2offresfile,s3files,rf] = findfiles(sourcedir)
    %% State 2 resonant imaging
    finds2img = questdlg('Do you have a state 2 image?', ...
        'Select state 2 image', ...
        'Yes','No','No');
    % Handle response
    switch finds2img
        case 'Yes'
            [s2file,s2path,~] = uigetfile('*.fits','Select the dataset',sourcedir,'MultiSelect','off');
        case 'No'
            s2file = '';
            s2path = '';
    end

    s2file = strcat(s2path,s2file);

    %% State 2 off-resonant image   
    finds2offresimg = questdlg('Do you have an off-resonant state 2 image?', ...
        'Select image', ...
        'Yes','No','No');
    % Handle response
    switch finds2offresimg
        case 'Yes'
            [s2offresfile,s2offrespath,~] = uigetfile('*.fits','Select the dataset',sourcedir,'MultiSelect','off');
        case 'No'
            s2offresfile = '';
            s2offrespath = '';
    end

    s2offresfile = strcat(s2offrespath,s2offresfile);

    %% State 3 images
    f = msgbox('Select all the state 3 images in the dataset');
    waitfor(f);
    [s3filenames,s3PathName,~] = uigetfile('*.fits','Select the dataset',sourcedir,'MultiSelect','on');
    s3files = strcat(s3PathName,s3filenames);
    
    for i=1:length(s3filenames)
        img_name = s3filenames{i};
        snipout = GetSnippetValues(img_name,{'RF23'});
        rf{i} = str2double(snipout.value{1});
    end
    
end







