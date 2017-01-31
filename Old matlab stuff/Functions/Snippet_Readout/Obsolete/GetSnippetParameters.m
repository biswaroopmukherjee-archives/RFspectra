function [parameter_value,match_timestamp,error_message] = ...
    GetSnippetValues(match_snippet_line,match_timestamp,error_message)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function matches a parameter line of the snippet file (located in
% the directory 'snippet_folder') to a given image (specified by 'img_name') within
% a time difference of a few seconds.
% The optional input varialble 'select_parameter' can be used to select a
% set of parameters for the output variable 'parameter_value'. If
% select_parameter is not specified all available parameters will be
% included in the output variable.
% If the program cannot find a match or finds more than one the 'error'
% will be true
%
% Input variables are of the following type: 
% string: filepath, img_name
% cell string: select_parameter
%
% Output variables are of the type:
% cell string: parameter_value
% string: snippet_timestamp, error_message
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Initialize required filepaths

user_folder = fileparts(fileparts(userpath));
dropbox_mit_BEC1 = '/Dropbox (MIT)/BEC1/';
snippet_folder = (strcat(user_folder,dropbox_mit_BEC1,'Image Data and Cicero Files/Data - Raw Images/Snippet_output/'));


%%% Initialize output variables

error_message = 'no error';
match_timestamp = '01-01-0000_00_00_00';

% if nargin<=2
%     select_parameter =
% end

parameter_value = NaN(size(select_parameter));


%%% Identify matching snippet file

% build snippet filename and convert time strings into serial numbers
timepart_img_name = img_name(1:19);

format_in = 'mm-dd-yyyy_HH_MM_ss';
img_name_serialnumber = datenum(timepart_img_name,format_in);

format_out = 'yyyy-mm-dd';
snippet_filename = [datestr(img_name_serialnumber,format_out) '.txt'];

% read the snippet file and match the timestamps of the image with the
% corresponding snippet line
[snippet_timestamp,parameter_string] = ReadSnippetFile(strcat(snippet_folder,snippet_filename));
snippet_timestamp_serialnumber = datenum(snippet_timestamp);

snippet_delay = snippet_timestamp_serialnumber - img_name_serialnumber;
snippet_lag = 3/(3600*24); % max lag of snippet timestamp after image in days
snippet_advance = -3/(3600*24); % max advance of snippet timestamp before image in days

match_vector_lag = snippet_delay <= snippet_lag;
match_vector_advance = snippet_delay >= snippet_advance;
match_vector = match_vector_lag.*match_vector_advance;


%%% Error handeling concerning matches (no match or more than one)

error_check_match = sum(match_vector);

if error_check_match == 0  %no match: try snippet file from next day
    
    %check next days 
    snippet_filename = [datestr(img_name_serialnumber+1,format_out) '.txt'];
    
    % read the snippet file and match the timestamps of the image with the
    % corresponding snippet line
    [timestamp,parameter_string] = ReadSnippetFile(strcat(snippet_folder,snippet_filename));
    timestamp_serialnumber = datenum(timestamp);

    snippet_delay = timestamp_serialnumber - img_name_serialnumber;

    match_vector_lag = snippet_delay <= snippet_lag;
    match_vector_advance = snippet_delay >= snippet_advance;
    match_vector = match_vector_lag.*match_vector_advance;
    
    error_check_match = sum(match_vector);
    
    if error_check_match == 1 % if there is a match with neighbouring days files
        error_message = 'Image matches with entry from neighbouring day''s snippet file.';
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
    matched_parameter_string = parameter_string(match_vector==1);
    
    % Regular expression to pick the values for the selected parameters.
    % Picks: after <parameter> and semicolon the stuff in parantheses(
    % optional +or- then [digits(1+).digits(1+)] or [digits(1+)]) before comma
    suffix_exp = ';([-+]?\d+\.\d+|\d+),';
    select_parameter_expression = strcat(select_parameter,suffix_exp);
    parameter_value = regexp(matched_parameter_string,select_parameter_expression,'tokens','once');
    
    % regexp outputs nested cells for the existing expression and empty cells
    % for the non-existing expression. The empty cells are converted into
    % 'NaN' nested cells and afterwards everything is flattened out 
    non_existing_parameter = cellfun('isempty',parameter_value);
    parameter_value(non_existing_parameter) = {{'NaN'}};
    parameter_value=[parameter_value{:}];
    
    if sum(non_existing_parameter)>0
        
        missing_parameters_string = strjoin(select_parameter(non_existing_parameter),', ');
        error_message=strcat('The following parameters do not exist: ',missing_parameters_string);
        
    end
    
end


end
