function [m] = mass_dz10(varargin)
if nargin
    if strcmp(varargin{1}, 'regions')
        load('mass_dz10_regions.mat', 'mass_t');
    else
        error('invalid input: only accepts "regions" or none');
    end
else
    load('mass_dz10.mat', 'mass_t');
end
% m = moving_average(mass_t, 15);
m = mass_t;
end
