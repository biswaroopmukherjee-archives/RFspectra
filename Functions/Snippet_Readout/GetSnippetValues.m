function [output,error_message] = ...
    GetSnippetValues(img_name,varargin)

% This function matches a parameter line of the snippet file (located in
% the directory 'snippet_folder') to a given image (specified by 
% 'img_name') within a time difference of a few seconds.
% The optional input varialble 'select_parameter' can be used to select a
% set of parameters for the output variable 'parameter_value'. If
% select_parameter is not specified all available parameters will be
% included in the output variable. Furthermore the optional input variable
% 'SnippetFolder' (needs to be called with 'SnippetFolder') can be 
% used to specify a non-default filepath for the snippet file.
%
% Input variables are of the following type: 
% string: img_name, snippet_filepath
% cell string: select_parameter
% 
% Example: GetSnippetValues(img_name,{'RF23'},'SnippetFolder','C:\')
%      or: GetSnippetValues(img_name,{'RF23'})
%      or: GetSnippetValues(img_name,'SnippetFolder','C:\')
%      or: GetSnippetValues(img_name)
%
% Output variables are of the type:
% cell string: parameter_value
% string: snippet_timestamp, error_message

%%% Initialize default snippet filepath  
user_folder = fileparts(fileparts(userpath));
dropbox_mit_BEC1 = '/Dropbox (MIT)/BEC1/';
default_snippet_folder = (fullfile(user_folder,dropbox_mit_BEC1,...
    'Image Data and Cicero Files/Data - Raw Images/Snippet_output/'));

%%% Initialize default selected parameters (for now RF23)
default_select_parameter = {'all parameters'};

% create input parser
p = inputParser;

addRequired(p,'ImageName',@ischar);
addOptional(p,'SelectParameter',default_select_parameter,@iscellstr)
addParameter(p,'SnippetFolder',default_snippet_folder,@ischar)

parse(p,img_name,varargin{:});

select_parameter = p.Results.SelectParameter;
snippet_folder = p.Results.SnippetFolder;

% get corresponding snippet string
[match_snippet_line,match_timestamp,error_message] = GetSnippetString(img_name,snippet_folder);

% initialze output parameter
parameter_value = cell(size(select_parameter));
parameter_value(:) = {'NaN'};

%%% Extracting the parameters for the matched snippet line
if strcmp(error_message,'no error') % match successful
    
        % if select_parameter is not specified all parameters are selected
        if strcmp(select_parameter,'all parameters')
            %%% Extracting the parameter names for the matched snippet line
            % Regular expression to pick the all strings in between a comma and a
            % semicolon. This misses the first parameter in the string!!
            expression = '\,(.*?)\;';
            all_parameters = regexp(match_snippet_line,expression,'tokens');
            % unnest cells
            select_parameter=[all_parameters{:}{:}];
            
            % get first parameter in string (not surrounded by , and ;)
            first_expression = '(.*?)\;';
            first_parameter = regexp(match_snippet_line,first_expression,'tokens','once');
            % unnest cell
            first_parameter = [first_parameter{1}];
            % combine parameters
            
            select_parameter=[first_parameter,select_parameter];
        end
    
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
    size_select_parameter = size(select_parameter);
    size_parameter_value = size(parameter_value);
    if size_select_parameter(1) ~= size_parameter_value(1)
       parameter_value = parameter_value';
    end
    
    if sum(non_existing_parameter)>0
        
        missing_parameters_string = strjoin(select_parameter(non_existing_parameter),', ');
        error_message=strcat('The following parameters do not exist: ',missing_parameters_string);
        
    end

    
end

% build a structure array out of select_parameter and parameter_value
  output.image = img_name;
  output.snippet = match_timestamp;
  output.parameter = select_parameter;
  output.value = parameter_value;

end
