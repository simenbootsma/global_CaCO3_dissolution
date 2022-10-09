function section35_diffusivity_uncertainty(data)
load('constants.mat');

%% Compute distributions
load('KD_distributions.mat', 'D_distribution', 'K_distribution'); % see compute_kd_distributions() below for details


%% Show distributions
figure;
subplot(2,1,1); hold on;
pbaspect([1 0.3072 0.3072]);
kd_facs = [0.5, 1.0, 1.5];
dcolors = [BLUES(2,:); PURPLES(2,:); REDS(2,:)];
for j = 1:length(kd_facs)
    plot(D_distribution{j,1}, D_distribution{j,2}, '-ok', 'MarkerFaceColor', dcolors(j,:), 'MarkerEdgeColor', dcolors(j,:), 'MarkerSize', 3);
end
xlabel('Isotropic diffusivity (m^2/s)', 'fontsize', 11, 'fontweight', 'bold');
ylabel('Volume fraction (-)', 'fontsize', 11, 'fontweight', 'bold');
xticks(-8:-1);
xticklabels(["10^{-8}", "10^{-7}", "10^{-6}", "10^{-5}", "10^{-4}", "10^{-3}","10^{-2}", "10^{-1}"])
ylim([0, 0.08]);
xlim([-8, -1]);
yticklabels(compose("%d%%", 0:2:8))
legend(compose("%.1f D", kd_facs), 'fontsize', 9);
grid on;
set(gcf, 'name', mfilename+" 1");


subplot(2,1,2); hold on;
pbaspect([1 0.3072 0.3072]);
kcolors = GREENS;
for j = 1:length(kd_facs)
    plot(K_distribution{j,1}, K_distribution{j,2}, '-ok', 'MarkerFaceColor', kcolors(j,:), 'MarkerEdgeColor', kcolors(j,:), 'MarkerSize', 3);
end
xlabel('Mesoscale diffusivity (m^2/s)', 'fontsize', 11, 'fontweight', 'bold');
ylabel('Volume fraction (-)', 'fontsize', 11, 'fontweight', 'bold');
xticks(-1:5);
xticklabels(["10^{-1}", "10^{0}", "10^{1}", "10^{2}", "10^{3}", "10^{4}", "10^{5}"])
legend(compose("%.1f K", kd_facs), 'fontsize', 9);
ylim([0, 0.08]);
xlim([-1, 5]);
yticklabels(compose("%d%%", 0:2:8))
grid on;
set(gcf, 'name', mfilename+" 2");



%% Show dissolution rate profiles
smth_fac = 1;
lw = 1.2;
strct = data.accum_alkP_kd;
colors = [[.32, .32, .32]; BLUES(1:3,:);PURPLES(1:3,:);REDS(1:3,:)];

[ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
set(gcf, 'name', mfilename+" 2");
if PER_UNIT_MASS
    xlim(ax1, [-.5 .5])
else
    xlim(ax0, [-200, 250])
    xlim(ax1, [-.25, .15])
end
xlim(ax2, xlim(ax1))
xticks(ax1, xticks(ax2))

x = strct(1,1).gamma;
y_without_bdy = [];
y = [];
for k = 1:3
    for d = 1:3
        y(3*d+k-3,:) = -0.5 * strct(k,d).tot;
        y_without_bdy(3*d+k-3,:) = -0.5 * (sum(strct(k,d).E(2:4,:), 1) + sum(strct(k,d).M, 1));
    end
end
y_smth = moving_average(y, smth_fac); % smoothen curves
y_without_bdy_smth = moving_average(y_without_bdy, smth_fac);

% Density plot
fill_y = [data.sulpis.gamma, fliplr(data.sulpis.gamma)];
fill_x = [data.sulpis.d_rate + data.sulpis.d_rate_std, fliplr(data.sulpis.d_rate - data.sulpis.d_rate_std)];
fill(ax0, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
plot(ax0, data.sulpis.d_rate, data.sulpis.gamma, "Color", [.3 .3 .3], "LineWidth", 1.2);

for k = 1:3
    for d = 1:3
        plot(ax0, y_smth(3*d+k-3,:), x, "Color", colors(3*d+k-2,:), "Linewidth", lw);
    end
end

% Depth plot
sulpis_depth_d_rate = data.sulpis.d_rate;
sulpis_depth_d_rate_std = data.sulpis.d_rate_std;
if PER_UNIT_MASS
    m = mass_for_densities(x);
    y = y .* m;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
    y_without_bdy = y_without_bdy .* m;
    sulpis_depth_d_rate =  sulpis_depth_d_rate .*  data.sulpis.mass;
    sulpis_depth_d_rate_std =  sulpis_depth_d_rate_std .*  data.sulpis.mass;
else
    y = y * 0.05;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
    y_without_bdy = y_without_bdy * 0.05;
    sulpis_depth_d_rate =  sulpis_depth_d_rate * 0.05;
    sulpis_depth_d_rate_std =  sulpis_depth_d_rate_std * 0.05;
end
load('depth_density_matrix_dz10.mat', 'Z', 'Gamma', 'A');
data.sulpis.gamma(isnan(data.sulpis.gamma)) = 0;
sulpis_depth_d_rate(isnan(sulpis_depth_d_rate)) = 0;
sulpis_depth_d_rate_std(isnan(sulpis_depth_d_rate_std)) = 0;
zt = sort(unique(Z(:)))';
B_sulpis = interp2(Z, Gamma, A, zt, data.sulpis.gamma', "linear", 0);
ys = sulpis_depth_d_rate * B_sulpis;
ys_min = (sulpis_depth_d_rate - sulpis_depth_d_rate_std) * B_sulpis;
ys_max = (sulpis_depth_d_rate + sulpis_depth_d_rate_std) * B_sulpis;

x(isnan(x)) = 0;
y(isnan(y)) = 0;
y_without_bdy(isnan(y_without_bdy)) = 0;
B = interp2(Z, Gamma, A, zt, x', "linear", 0);
y_without_bdy = y_without_bdy * B;


if PER_UNIT_MASS
    m = mass_dz10();
%     y = y ./ m; % now divide by mass of each depth layer to return to units of umol/kg/yr
    ys = ys ./ m;
    ys_min = ys_min ./ m;
    ys_max = ys_max ./ m;
    y_without_bdy = y_without_bdy ./ m;
else
%     y = y ./ 10; % divide by dz
    ys = ys ./ 10;
    ys_min = ys_min ./ 10;
    ys_max = ys_max ./ 10;
    y_without_bdy = y_without_bdy ./ 10;
end

fill_y = [zt, fliplr(zt)] * 1e-3;
fill_x = [ys_max, fliplr(ys_min)];

fill(ax1, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
plot(ax1, ys, zt*1e-3, "Color", [.3 .3 .3], "LineWidth", 1.2);
fill(ax2, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
plot(ax2, ys, zt*1e-3, "Color", [.3 .3 .3], "LineWidth", 1.2);

for k = 1:3
    for d = 1:3
        plot(ax1, y_without_bdy(3*d+k-3,:), zt*1e-3, "Color", colors(3*d+k-2,:), "LineWidth", lw);
        plot(ax2, y_without_bdy(3*d+k-3,:), zt*1e-3, "Color", colors(3*d+k-2,:), "LineWidth", lw);
    end
end

axes('position', [.55, .33, .15, .15]);
[im, map, alpha] = imread('diff_unc_inset.png');
f = imshow(im, map);
legend(ax2, ["", "", "Alk^*"], "fontsize", 9, "position", [0.575892247468979,0.271698310971386,0.129166666666667,0.053968253968254]);
set(f, 'alphadata', alpha);
set(gcf, 'name', mfilename+" 3");

end




function compute_kd_distributions()
    load("Casimir_D_Tgrid.mat", "Casimir_D_Tgrid");  % D estimates from Casimir DeLavergne
    load("WOA18_K.mat", "WOA18_K"); % K estimates from World Ocean Atlas
    load("WOA18_Volume.mat", 'vol');
    D = mean(Casimir_D_Tgrid.Dtot, 4, 'omitnan');
    K = WOA18_K.K_WOA(:,:,2:end);
    vol(isnan(D+K)) = nan;
    tot_vol = sum(vol(:), 'omitnan');
    
    kdfacs = [0.5, 1.0, 1.5];
    D_distribution = cell(length(kdfacs),2);
    K_distribution = cell(length(kdfacs),2);
    for j = 1:length(kdfacs)
        [~, edges_d, idx_d] = histcounts(log10(kdfacs(j) * D), 'binwidth', 0.1);
        x_d = edges_d(1:end-1) + diff(edges_d)/2;
        y_d = nan(1, length(x_d));
        for i = 1:length(x_d)
            y_d(i) = sum(vol(idx_d==i), 'omitnan') / tot_vol;
        end
        D_distribution(j,:) = {x_d, y_d};
        
        [~, edges_k, idx_k] = histcounts(log10(kdfacs(j) * K), 'binwidth', 0.1);
        x_k = edges_k(1:end-1) + diff(edges_k)/2;
        y_k = nan(1, length(x_k));
        for i = 1:length(x_k)
            y_k(i) = sum(vol(idx_k==i), 'omitnan') / tot_vol;
        end
        K_distribution(j,:) = {x_k, y_k};
    end
    save('KD_distributions.mat', 'D_distribution', 'K_distribution');
end
