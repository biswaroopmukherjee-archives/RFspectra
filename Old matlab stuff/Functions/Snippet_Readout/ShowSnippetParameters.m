function [all_parameters,varied_parameters,error_message] = ...
    ShowSnippetParameters(start_img_name,end_image_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function shows for a single image or a series of images all
% parameter names and the parameters that have been changed inbetween
% different images of the series
%
% Input variables are of the following type: 
% string: start_img_name, end_image_name (optional)
%
% Output variables are of the type:
% cell string: all_parameters, varied_parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Initialize required filepaths
year_img = start_img_name(7:10);
month_img = strcat(year_img,'-',start_img_name(1:2));
day_img = strcat(month_img,'-',start_img_name(4:5));
user_folder = fileparts(fileparts(userpath));
dropbox_mit_BEC1 = '/Dropbox (MIT)/BEC1/';
img_subfolder = strcat('Image Data and Cicero Files/Data - Raw Images/',year_img,...
    '/',month_img,'/',day_img,'/');
img_folder = fullfile(user_folder,dropbox_mit_BEC1,img_subfolder);


%%% get the snippet strings for each image
[match_snippet_line,match_timestamp,error_message] = GetSnippetString(start_img_name);

%%% Extracting the parameter names for the matched snippet line
% Regular expression to pick the all strings in between a comma and a
% semicolon
expression = '\,(.*?)\;';
all_parameters = regexp(match_snippet_line,expression,'tokens');
% unnest cells
all_parameters=[all_parameters{:}{:}];

end
