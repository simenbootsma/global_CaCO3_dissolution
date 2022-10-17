function VOL_JDIST = volume_jdist(vol,var1,var2,var1_range,var2_range)

VOL_JDIST = zeros(length(var1_range)-1,length(var2_range)-1);

%% ADDED BY SJOERD ON 24-03-2018
% REMOVE NANS for faster binning
idx = isnan(vol);
vol(idx)  = [];
var1(idx) = [];
var2(idx) = [];

for onedx = 1:length(var1_range)-1
  %find all fluxes with value var1 and bin
  %them according to their var2 value
  [~,BIN] = histc(var2(var1>=var1_range(onedx) &...
             var1<var1_range(onedx+1)),var2_range);
  BIN(BIN==0)=size(var2_range,2)+1; %because it doesn't like zero
  %Sum all of these fluxes 
  if isfinite(BIN)
    totfreq = accumarray(BIN,...
    double(vol(var1>=var1_range(onedx) &...
    var1<var1_range(onedx+1))),...
    [length(var2_range)+1 1]);
    VOL_JDIST(onedx,:) = totfreq(1:end-2);
  end
end

% - Imagine input: 
% - volume_jdist(H,S,T,Sedge,Tedge)
% - Then we have:
% - histc(T(S>=Sedge(i) & S<Sedge(i+1)),Tedge);
% - First the temperatures in a certain salinity range are selected.
% - Then these temperatures are catagorized in bins that span the whole temperature domain.
% - Bin thus gives the indexes of all temperatures (T) that are in a
% - certain temperature and salinity range.
% - Or: All temperatures within a certain salinity range are binned in
% - temperature bins. BINS gives the indexes of vector
% - T(S>=Srange(i) & S<Srange(i+1)) that are within a certain T range for
% - each T range.
% - accumarray(BIN,H(S>=Srange(i) & S<Srange(i+1)),[Length(Trange)+1 1]);
% - All heatfluxes within a certain salinity range, that are in a certain 
% - temperature range are summed. Such that, for each temperature range
% - (and salinity grid) one has the total heat flux.
% - Iterating this over all salinity ranges, this means that for each dTdS
% - increament, the total heat flux in that area is summed.