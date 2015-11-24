function plotspec(spec,rf)
%% PLOTSPEC plots an arbitrary spectrum image given a cell array of rf values
imagesc(spec);
set(gca,'XTickLabel',num2str(cell2mat(rf(2:2:end))));
end