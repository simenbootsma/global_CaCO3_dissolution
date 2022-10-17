function binned = binning_1D(values_to_bin,label_for_binning,bin_edges)

% HERE WE BIN "values_to_bin" ACCORDING TO THE VALUES IN
% "label_for_binning", FOR BINS OF THE SIZE "bin_edges".

values_to_bin(isnan(values_to_bin))=0; % Remove NaN-values.
[~,idx]     = histc(label_for_binning(:),bin_edges); % Obtain index determining which bin it belongs.
idx0        = idx~=0; % idx = 0 where label_for_binning=NaN, so this gives the indeces where label_for_binning = not 0.
idx         = idx(idx0); % Select only where idx>0;
binned      = accumarray(idx,values_to_bin(idx0),[length(bin_edges)-1 1]); % Do the binning, an maintain bin-length vector.
binned      = binned';

% function binned = binning_1D(values_to_bin,label_for_binning,bin_edges)
% 
% % HERE WE BIN "values_to_bin" ACCORDING TO THE VALUES IN
% % "label_for_binning", FOR BINS OF THE SIZE "bin_edges".
% 
% values_to_bin(isnan(values_to_bin))=0; % Remove NaN-values.
% [~,~,idx]     = histcounts(label_for_binning(:),bin_edges); % Obtain index determining which bin it belongs.
% idx0        = idx~=0; % idx = 0 where label_for_binning=NaN, so this gives the indeces where label_for_binning = not 0.
% idx         = idx(idx0); % Select only where idx>0;
% binned      = accumarray(idx,values_to_bin(idx0),[length(bin_edges)-1 1]); % Do the binning, an maintain bin-length vector.
% binned      = binned';

