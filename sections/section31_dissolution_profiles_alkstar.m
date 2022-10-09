function section31_dissolution_profiles_alkstar(data)
load('constants.mat');

%% Total
[ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
if PER_UNIT_MASS
    xlim(ax0, [-1 1]);
    xlim(ax1, [-.2 .5]);
else
    xlim(ax0, [-100, 200]);
    xlim(ax1, [-0.05, 0.12]);
end
xlim(ax2, xlim(ax1));
xticks(ax1, xticks(ax2));

% density plot
fill_y = [data.sulpis.gamma, fliplr(data.sulpis.gamma)];
fill_x = [data.sulpis.d_rate + data.sulpis.d_rate_std, fliplr(data.sulpis.d_rate - data.sulpis.d_rate_std)];

fill(ax0, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
plot(ax0, data.sulpis.d_rate, data.sulpis.gamma, "Color", [.3 .3 .3], "LineWidth", 1.2);

% depth plot
if PER_UNIT_MASS
     data.sulpis.d_rate =  data.sulpis.d_rate .*  data.sulpis.mass;
     data.sulpis.d_rate_std =  data.sulpis.d_rate_std .*  data.sulpis.mass;
else
    % multiply by bin width
     data.sulpis.d_rate =  data.sulpis.d_rate .* 0.05;
     data.sulpis.d_rate_std =  data.sulpis.d_rate_std .* 0.05;
end
load('depth_density_matrix_dz10.mat', 'Z', 'Gamma', 'A');
data.sulpis.gamma(isnan(data.sulpis.gamma)) = 0;
data.sulpis.d_rate(isnan(data.sulpis.d_rate)) = 0;
zt = sort(unique(Z(:)))';
B = interp2(Z, Gamma, A, zt, data.sulpis.gamma', "linear", 0);
y = data.sulpis.d_rate * B;
y_min = (data.sulpis.d_rate - data.sulpis.d_rate_std) * B;
y_max = (data.sulpis.d_rate + data.sulpis.d_rate_std) * B;
if PER_UNIT_MASS
    m = mass_dz10();
    y = y ./ m;
    y_min = y_min ./ m;
    y_max = y_max ./ m;
else
    % divide by bin width
    y = y / 10;
    y_min = y_min / 10;
    y_max = y_max / 10;
end
fill_y = [zt, fliplr(zt)] * 1e-3;
fill_x = [y_max, fliplr(y_min)];

fill(ax1, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
plot(ax1, y, zt*1e-3, "Color", [.3 .3 .3], "LineWidth", 1.2);
fill(ax2, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
plot(ax2, y, zt*1e-3, "Color", [.3 .3 .3], "LineWidth", 1.2);

set(gcf, 'name', mfilename+" 1");





%% Per region
[ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
if PER_UNIT_MASS
    xlim(ax0, [-1 1]);
    xlim(ax1, [-.2 .5]);
else
    xlim(ax0, [-10, 20]);
    xlim(ax1, [-.005 .012]);
end
xlim(ax2, xlim(ax1));
xticks(ax1, xticks(ax2));

% density plot
for i = 1:size(data.sulpis.d_rate_biomes, 1)
    plot(ax0, data.sulpis.d_rate_biomes(i,:), data.sulpis.gamma, "Color", BIOME_COLORS(i,:), "LineWidth", 1.5);
end

% depth plot
if PER_UNIT_MASS
    data.sulpis.d_rate_biomes =  data.sulpis.d_rate_biomes .*  data.sulpis.biomes_mass;
else
    data.sulpis.d_rate_biomes =  data.sulpis.d_rate_biomes .* 0.05;
end
load('depth_density_matrix_dz10.mat', 'Z', 'Gamma', 'A');
data.sulpis.gamma(isnan(data.sulpis.gamma)) = 0;
data.sulpis.d_rate_biomes(isnan(data.sulpis.d_rate_biomes)) = 0;
zt = sort(unique(Z(:)))';
B = interp2(Z, Gamma, A, zt, data.sulpis.gamma', "linear", 0);
y = data.sulpis.d_rate_biomes * B;
if PER_UNIT_MASS
    m = mass_dz10('regions');
    m(m < 1e16) = nan;
    y = y ./ m;
else
    y = y / 10;
end

for i = 1:size(data.sulpis.d_rate_biomes, 1)
    plot(ax1, y(i,:), zt*1e-3, "Color", BIOME_COLORS(i,:), "LineWidth", 1.5);
    plot(ax2, y(i,:), zt*1e-3, "Color", BIOME_COLORS(i,:), "LineWidth", 1.5);
end

colormap(BIOME_COLORS);
cb = colorbar(ax0, 'location', 'NorthOutside');

cb.Ticks = .05:.1:1;
cb.TickLabels = 1:10;
cb.Label.String = "Region";
cb.Label.FontSize = 11;
cb.Label.FontWeight = "bold";
cb.Position = [0.16577380952381,0.88047619047621,0.67113095238095,0.030793650793651];
set(gcf, 'name', mfilename+" 2");

end