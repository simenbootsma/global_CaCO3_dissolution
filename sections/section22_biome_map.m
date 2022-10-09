function section22_biome_map(~)
load('constants.mat');
load("GLODAPv2_Biomes2D.mat", "biomes2d", "biomes_extended2d");
load("WOA18_Coordinates.mat", "x", "y");


% compute map
x2d = repmat(x', [1, length(y)]);
y2d = repmat(y, [length(x), 1]);

biomes2d(isnan(biomes_extended2d)) = 0;
biomes2d = biomes2d';
biomes2d = [biomes2d(:, 181:end), biomes2d(:, 1:180)];

biomes_extended2d(isnan(biomes_extended2d)) = 0;
biomes_extended2d = biomes_extended2d';
biomes_extended2d = [biomes_extended2d(:, 181:end), biomes_extended2d(:, 1:180)];

I0 = isnan(biomes2d) & (mod((x2d-y2d)',10) < 2); % stripes
I1 = isnan(biomes2d) & (biomes2d([end, 1:end-1], :) > 0); % up
I2 = isnan(biomes2d) & (biomes2d([2:end, 1], :) > 0); % down
I3 = isnan(biomes2d) & (biomes2d(:, [end, 1:end-1]) > 0); % left
I4 = isnan(biomes2d) & (biomes2d(:, [2:end, 1]) > 0); % right
biomes_extended2d(logical(I0+I1+I2+I3+I4)) = -1;
extended_biomes_map = biomes_extended2d;

% show figure
fig = figure; hold on;
imagesc([-180, 180], [-90, 90], extended_biomes_map)
xlabel('Longitude', 'fontsize', 11, 'fontweight', 'bold');
ylabel('Latitude', 'fontsize', 11, 'fontweight', 'bold');
pbaspect([2 1 1])
xlim([-180, 180]);
ylim([-90, 90]);
xticks([-180, -120, -60, 0, 60, 120, 180]);
yticks([-90, -60, -30, 0, 30, 60, 90]);
xticklabels(["180"+sprintf("%c", char(176))+"W", "120"+sprintf("%c", char(176))+"W", "60"+sprintf("%c", char(176))+"W", "0"+sprintf("%c", char(176)), ...
    "60"+sprintf("%c", char(176))+"E", "120"+sprintf("%c", char(176))+"E","180"+sprintf("%c", char(176))+"E"]);
yticklabels(["90"+sprintf("%c", char(176))+"S", "60"+sprintf("%c", char(176))+"S", "30"+sprintf("%c", char(176))+"S", "0"+sprintf("%c", char(176)), ...
    "30"+sprintf("%c", char(176))+"N", "60"+sprintf("%c", char(176))+"N","90"+sprintf("%c", char(176))+"N"]);

colormap([[.7,.7,.7;0, 0, 0]; BIOME_COLORS]);
caxis([-1.5, 10.5]);

cb = colorbar;
labels = ["1: Subpolar North Pacific", "2: Subtropical North Pacific", "3: Equatorial Pacific", "4: Subtropical South Pacific", "5: Subpolar North Atlantic", "6: Subtropical North Atlantic", "7: Equatorial Atlantic", "8: Subtropical South Atlantic", "9: Indian Ocean", "10: Southern Ocean"];
cb.Ticks = 1:10;
cb.TickLabels = labels;

fig.Position(3) = 1010;
cb.Limits = [.5, 10.5];
fig.CurrentAxes.Position(1) = 0.08;
set(gcf, 'name', mfilename);

end