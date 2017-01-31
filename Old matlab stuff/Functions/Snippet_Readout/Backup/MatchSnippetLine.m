function [match_vector] = MatchSnippetLine(img_name,snippet_timestamp,lag,advance)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds matching snippet lines within a certain 'lag' and
% 'advance' interval around the timestamp of the image.
% The ouput variable is 'match_vecor' which has entries of 1 when there is a
% match and zero for no match.
% 
% Input variables are of the following type: 
% string: img_name
% double: lag, advance (in seconds)
% cellstring : snippet_timestamp (all timestamps from a snippet file)
%
% Output variables are of the type:
% string: match_vec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% convert snippet and image timestamp into unique serialnumbers
snippet_timestamp_serialnumber = datenum(snippet_timestamp);
timepart_img_name = img_name(1:19);
format_in = 'mm-dd-yyyy_HH_MM_ss';
img_name_serialnumber = datenum(timepart_img_name,format_in);

% compare the serialnumbers to get the match
snippet_delay = snippet_timestamp_serialnumber - img_name_serialnumber;
snippet_lag = lag/(3600*24); % max lag of snippet timestamp after image in days
snippet_advance = -advance/(3600*24); % max advance of snippet timestamp before image in days

match_vector_lag = snippet_delay <= snippet_lag;
match_vector_advance = snippet_delay >= snippet_advance;
match_vector = match_vector_lag.*match_vector_advance;

end