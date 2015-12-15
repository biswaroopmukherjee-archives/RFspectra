function [match_snippet_line,match_timestamp,error_message] = GetSnippetString(img_name)

% This function matches a line (string) of the snippet file (located in
% the directory 'snippet_folder') to a given image (specified by 
%'img_name') within a time difference of a few seconds. The default 
% filepath for the snippet file is '.../Dropbox (MIT)/BEC1/Image Data and \
% Cicero Files/Data - RawImages/Snippet_output/'. Optionally another path
% can be specified with the input variable 'snippet_filepath'.
%
% Input variables are of the following type: 
% string: img_name, snippet_filepath
%
% Output variables are of the type:
% string: match_snippet_line, snippet_timestamp, error_message


%%% Initialize required filepaths   
user_folder = fileparts(fileparts(userpath));
dropbox_mit_BEC1 = '/Dropbox (MIT)/BEC1/';
snippet_folder = (fullfile(user_folder,dropbox_mit_BEC1,'Image Data and Cicero Files/Data - Raw Images/Snippet_output/'));

%%% Initialize output variables

error_message = 'no error';
match_timestamp = '01-01-0000_00_00_00';
match_snippet_line = 'empty';

%%% Set tolerances for snippet match in seconds
lag = 3;
advance = 3;

%%% Identify matching snippet file

% build snippet filename and convert time strings into serial numbers
timepart_img_name = img_name(1:19);

format_in = 'mm-dd-yyyy_HH_MM_ss';
img_name_serialnumber = datenum(timepart_img_name,format_in);

format_out = 'yyyy-mm-dd';
snippet_filename = [datestr(img_name_serialnumber,format_out) '.txt'];

% read the snippet file and match the timestamps of the image with the
% corresponding snippet line
[snippet_timestamp,snippet_lines] = ReadSnippetFile(fullfile(snippet_folder,snippet_filename));
[match_vector] = MatchSnippetLine(img_name,snippet_timestamp,lag,advance);


%%% Error handeling concerning matches (no match or more than one)

error_check_match = sum(match_vector);

if error_check_match == 0  %no match: try snippet file from next day
    
    for i=-1:2:1 % for loop to check neighbouring day's snippet files
    
    %check next days
    snippet_filename = [datestr(img_name_serialnumber + i,format_out) '.txt'];
    
    % check that the snippet file for the next days exists
    if exist(strcat(snippet_folder,snippet_filename), 'file') == 2
    
    % read the snippet file and match the timestamps of the image with the
    % corresponding snippet line
    [snippet_timestamp,snippet_lines] = ReadSnippetFile(strcat(snippet_folder,snippet_filename));
    [match_vector] = MatchSnippetLine(img_name,snippet_timestamp,lag,advance);
    
    error_check_match = sum(match_vector);
    
    if error_check_match == 1; % terminate for loop when finding unique match
        break
    end
    
    end
    
    end
    
    
    % Setting error messages
    if error_check_match == 1 % if there is a match with neighbouring days files
        error_message = 'no error';
    elseif error_check_match == 0 % if there is no match with neighbouring days files
        error_message = 'No match found in same day and neighbouring day''s snippet files.';
    elseif error_check_match > 1 % if there is a match with neighbouring days files but the match is not unique
        error_message = 'Found match in neighbouring day''s snippet file but the match is not unique.';
    end
    
    
elseif error_check_match > 1 % match not unique
    error_message = 'The match is not unique.';

end

%%% Extracting the parameters for the matched snippet line
if error_check_match == 1 % match successful

    % output matched snippet timestamp
    match_timestamp = snippet_timestamp(match_vector==1);
    
    % parameter string from the matching snippet line
    match_snippet_line = snippet_lines(match_vector==1);

end

end