function [m] = mass_for_densities(nd)
load("D:\NIOZ\Tracer data\WOA18_Density.mat", 'rho');
load("D:\NIOZ\Tracer data\WOA18_Volume.mat", 'vol');
load("D:\NIOZ\Tracer data\WOA18_NeutralDensity.mat", 'gamma');
addpath("D:\NIOZ\WMT\"); % binning_1D()
mass = vol .* (1000 + rho);
nd_edges = [nd - (nd(2)-nd(1))/2, nd(end) + (nd(2)-nd(1))/2]; % assumes density is equally spaced
m = binning_1D(mass, gamma, nd_edges) / size(mass, 4);
end
