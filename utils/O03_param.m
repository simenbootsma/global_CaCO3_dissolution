function [y] = O03_param(zq)
% load('D:\NIOZ\Tracer data\WOA18_ChlorophyllA.mat', 'Chla');
% Chla4d = permute(repmat(Chla,[1,1,1,zl]),[1,2,4,3]); % 4D Chla, 
z = zq(1:end-1) + diff(zq)/2;
z = [zq(1) - (zq(2)-zq(1))/2, z, zq(end) + (zq(end)-zq(end-1))/2];

%% Ohlman 2003:
mean_Chla = 0.16; % mg/m3, from Groeskamp et al. 2018
A1 	= 0.571 + 0.025 .*  log(0.149 .* mean_Chla);
A2	= 0.223 + 0.010 .*  log(2.329 .* mean_Chla);
B1	= 0.015 + 0.176 .* sqrt(0.462 .* mean_Chla);
B2	= 0.688 + 0.060 .*  log(0.125 .* mean_Chla);
fac_O03            = A1 .* exp(-B1 .* z) + A2 .* exp(-B2 .* z);
accum = diff(fac_O03); % heat accumulated between depth layers 
y = accum / sum(accum, 'omitnan');  % normalize
end