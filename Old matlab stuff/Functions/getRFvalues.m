function rf = getRFvalues(images)

%% Get filenames
for i=1:length(images)
    fullpath = images{i};
    f1 = findstr(fullpath,'/');
    f2 = findstr(fullpath,'\');
    f=[f1 f2];
    if isempty(f)
        startpoint = 1;
    else
        startpoint = 1+max(f);
    end
    filenames{i} = fullpath(startpoint:end-5);
end

for i=1:length(filenames)
        img_name = filenames{i};
        snipout = GetSnippetValues(img_name,{'RF23'});
        rf{i} = str2double(snipout.value{1});
end