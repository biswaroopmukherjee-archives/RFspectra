function [parameter_value,match_timestamp,error_message] = ...
    GetSnippetValues_test(img_name,varargin)

% This function matches a parameter line of the snippet file (located in
% the directory 'snippet_folder') to a given image (specified by 
% 'img_name') within a time difference of a few seconds.
% The optional input varialble 'select_parameter' can be used to select a
% set of parameters for the output variable 'parameter_value'. If
% select_parameter is not specified all available parameters will be
% included in the output variable. Furthermore the optional input variable
% 'SnippetFilepath' (needs to be called with 'snippet_filepath') can be 
% used to specify a non-default filepath for the snippet file.
%
% Input variables are of the following type: 
% string: img_name, snippet_filepath
% cell string: select_parameter
%
% Output variables are of the type:
% cell string: parameter_value
% string: snippet_timestamp, error_message

% Default path for snippet output file
user_folder = fileparts(fileparts(userpath));
dropbox_mit_BEC1 = '/Dropbox (MIT)/BEC1/';
snippet_folder = (fullfile(user_folder,dropbox_mit_BEC1,...
    'Image Data and Cicero Files/Data - Raw Images/Snippet_output/'));

% use input parser to address variable input arguments
% valid parameters, and their default values
    parameter_names = {'SelectParameters','SnippetFilepath'};
    default_values  = {'RF23',snippet_folder};

% parse function arguments
    [select_parameter,snippet_filepath] = ...
        internal.stats.parseArgs(parameter_names, default_values, varargin{:});
select_parameter
snippet_filepath
% use the processed values: clr, lw, ls, txt
% corresponding to the specified parameters

% get corresponding snippet string
[match_snippet_line,match_timestamp,error_message] = GetSnippetString(img_name,snippet_filepath);
match_snippet_line
% initialze output parameter
parameter_value = cell(size(select_parameter));
parameter_value(:) = {'NaN'};

%%% Extracting the parameters for the matched snippet line
if strcmp(error_message,'no error') % match successful
    
    % Regular expression to pick the values for the selected parameters.
    % Picks: after <parameter> and semicolon the stuff in parantheses(
    % optional +or- then [digits(1+).digits(1+)] or [digits(1+)]) before comma
    %suffix_exp = ';([-+]?\d+\.\d+|\d+),';
    suffix_exp='\;(.*?)\,';
    select_parameter_expression = strcat(select_parameter,suffix_exp);
    parameter_value = regexp(match_snippet_line,select_parameter_expression,'tokens','once');
    
    % regexp outputs nested cells for the existing expression and empty cells
    % for the non-existing expression. The empty cells are converted into
    % 'NaN' nested cells and afterwards everything is flattened out 
    non_existing_parameter = cellfun('isempty',parameter_value);
    parameter_value(non_existing_parameter) = {{'NaN'}};
    parameter_value=[parameter_value{:}];
    
    % transpose parameter_values if necessary to match with
    % select_parameters
    size_select_parameters = size(select_parameters);
    if size_select_parameters(2) > size_select_parameters(1)
       parameter_value = parameter_value';
    end
    
    if sum(non_existing_parameter)>0
        
        missing_parameters_string = strjoin(select_parameter(non_existing_parameter),', ');
        error_message=strcat('The following parameters do not exist: ',missing_parameters_string);
        
    end
    
end


end
