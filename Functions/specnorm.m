function specout = specnorm(spec)
%% SPECNORM normalizes the spectra
    s=size(spec);
    spec_sum = repmat(sum(spec,2),[1,s(2)]);
    specout = spec./spec_sum;
end
