function [m] = mass_for_densities_ext_regions(nd)
load("D:\NIOZ\Tracer data\WOA18_Density.mat", 'rho');
load("D:\NIOZ\Tracer data\WOA18_Volume.mat", 'vol');
load("D:\NIOZ\Tracer data\WOA18_NeutralDensity.mat", 'gamma');
load("D:\NIOZ\Tracer data\GLODAPv2_Biomes4D_extended.mat", 'biomes')
addpath("D:\NIOZ\WMT\"); % binning_1D()
mass = vol .* (1000 + rho);
nd_edges = [nd - (nd(2)-nd(1))/2, nd(end) + (nd(2)-nd(1))/2]; % assumes density is equally spaced

m = nan(10, length(nd_edges)-1);
for i = 1:10
    I = biomes == i;
    m(i,:) = binning_1D(mass(I), gamma(I), nd_edges) / size(mass, 4);
end
end
