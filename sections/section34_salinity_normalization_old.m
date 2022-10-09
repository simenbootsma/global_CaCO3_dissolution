function section34_salinity_normalization(data)
    load('constants.mat');
    load('D:\NIOZ\Tracer data\GLODAPv2_Biomes4D_extended.mat', 'biomes');
    load('D:\NIOZ\Tracer data\CANYONB_Alkalinity.mat', 'alk');
    load('D:\NIOZ\Tracer data\WOA18_AbsoluteSalinity.mat', 'SA');
    load('D:\NIOZ\Tracer data\WOA18_Nitrate.mat', 'Ni');
    load('D:\NIOZ\Tracer data\WOA18_Volume.mat', 'vol');
    load('D:\NIOZ\Tracer data\WOA18_Density.mat', 'rho');
    load('D:\NIOZ\Tracer data\WOA18_Coordinates4D.mat', 'z');
    mass = vol .* (1000 + rho);
    
    %% Fits per region
    surf_height = 20; % [m]
    n_biomes = length(unique(biomes(~isnan(biomes))));
    slopes = nan(1, n_biomes);
    slopes_std = nan(1, n_biomes);
    intercepts = nan(1, n_biomes);
    intercepts_std = nan(1, n_biomes);
    
    figure;
    names = ["Subpolar"+newline+"North Pacific", "Subtropical"+newline+"North Pacific", "Equatorial"+newline+"Pacific", "Subtropical"+newline+"South Pacific", "Subpolar"+newline+"North Atlantic", "Subtropical"+newline+"North Atlantic", "Equatorial"+newline+"Atlantic", "Subtropical"+newline+"South Atlantic", "Indian"+newline+"Ocean", "Southern"+newline+"Ocean"];
    fits = struct('name', [], 'salinity', [], 'alkp', [], ...
        'intercept', [], 'intercept_std', [], 'slope', [], 'slope_std', [], ...
        'global_intercept', [], 'global_intercept_std', [], 'global_slope', [], ...
        'global_slope_std', [], 'global_average', [], 'global_average_std', []);

    for j = 1:n_biomes
        I = z <= surf_height & biomes == j & ~isnan(alk);
        [b, sb] = lscov([ones(size(SA(I))), SA(I)], alk(I)+1.26*Ni(I), mass(I));
        intercepts(j) = b(1);
        intercepts_std(j) = sb(1);
        slopes(j) = b(2);
        slopes_std(j) = sb(2);
    
        subtightplot(2,5,j,0.02, [.1, .01], [.07 .01]);
        hold on;
        box on;
        scatter(SA(I), alk(I)+1.26*Ni(I), 7, BIOME_COLORS(j,:), 'filled', 'MarkerFaceAlpha', .1);
        plot([0, 45], b(2)*[0, 45] + b(1), '-k', 'linewidth', .8);
        if j > 5
            xlabel('S_A (g/kg)', 'fontsize', 13, 'fontweight', 'bold');
        else
            xticklabels([]);
        end
        if j == 1 || j == 6
            ylabel('A_P ({\mu}mol kg^{-1})', 'fontsize', 13, 'fontweight', 'bold');
        else
            yticklabels([]);
        end
    
        grid on;
        xlim([0, 45]);
        ylim([1000, 2700])
        xticks(0:15:45);
        yticks(1000:500:3000);
        set(gca, 'color', [.7 .7 .7]);
        set(gca, 'GridColor', [1 1 1]);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'FontSize', 11);
        text(2, 2600, "Region "+string(j), 'FontSize', 14, 'Color', BIOME_COLORS(j,:), 'FontWeight','bold')
        text(2, 2420, "$\alpha$ = "+sprintf("%.1f", b(2)), 'FontSize', 11, 'Color', 'k', 'interpreter', 'latex')
        text(2, 2300, "$\beta$ = "+sprintf("%.1f", b(1)), 'FontSize', 11, 'Color', 'k', 'interpreter', 'latex')
    
        fits(j).name = names(j);
        fits(j).salinity = SA(I);
        fits(j).alkp = alk(I) + 1.26 * Ni(I);
        fits(j).intercept = intercepts(j);
        fits(j).intercept_std = intercepts_std(j);
        fits(j).slope = slopes(j);
        fits(j).slope_std = slopes_std(j);
    end
    set(gcf, 'position', [1,41,1280,607.3333333333333]);
    set(gcf, 'name', mfilename+" 1");
    if SAVE_FIGS
        filename = "sec34_regional_fits";
        folder = "../auto_saved/";
        if PER_UNIT_MASS, folder = folder + "per_unit_mass/"; else folder = folder + "total_amount/"; end
        saveas(gcf, folder + filename + ".png");
        print(gcf, folder + filename + ".eps", '-depsc');
    end
    



    %% All fits combined
    % Global averages
    I = z <= surf_height & ~isnan(alk + Ni + SA);
    avg_alkP = sum(mass(I) .* (alk(I) + 1.26*Ni(I))) / sum(mass(I));
    avg_SA = sum(mass(I) .* SA(I)) / sum(mass(I));
    disp("Average AlkP: "+sprintf("%.1f", avg_alkP)+" umol/kg");
    disp("Average SA:   "+sprintf("%.2f", avg_SA)+" g/kg");
    
    % Global fits
    [glob_b, glob_sb] = lscov([ones(size(SA(I))), SA(I)], alk(I)+1.26*Ni(I), mass(I));
    disp("Global AlkP/Sal-slope:     "+string(glob_b(2))+" +/- "+string(glob_sb(2)));
    disp("Global AlkP/Sal-intercept: "+string(glob_b(1))+" +/- "+string(glob_sb(1)));
    
    [glob_b2, glob_sb2] = lscov(SA(I), alk(I)+1.26*Ni(I), mass(I));
    disp("Global AlkP/Sal-avg:     "+string(glob_b2(1))+" +/- "+string(glob_sb2(1)));

    fits(1).global_intercept = glob_b(1);
    fits(1).global_intercept_std = glob_sb(1);
    fits(1).global_slope = glob_b(2);
    fits(1).global_slope_std = glob_sb(2);
    fits(1).global_average = glob_b2(1);
    fits(1).global_average_std = glob_sb2(1);

    save('data\salinity_normalization_data.mat', 'fits')

    
    % Fig
    figure; hold on;
    I = find(z <= surf_height);
    scatter(SA(I),alk(I)+1.26*Ni(I),5,[.7 .7 .7],'filled','markerfacealpha',1);
    sal_x = 0:45;
    alk_y = slopes'*sal_x + intercepts';
    for i = 1:n_biomes
        plot(sal_x, alk_y(i,:),'color',BIOME_COLORS(i,:), 'linewidth',1);
    end
    h(1)=plot(sal_x, glob_b2(1)*sal_x, 'color','k', 'linewidth',1.5);
    h(2)=plot(sal_x, glob_b(2)*sal_x + glob_b(1), '--', 'color', 'k', 'linewidth',1.5);
    colormap(BIOME_COLORS(1:n_biomes,:));
    cb = colorbar;
    cb.Label.String = 'Region';
    cb.Label.FontSize = 11;
    cb.Label.FontWeight = 'bold';
    caxis([0.5, 10.5 ]);
    xlabel('S_A (g/kg)', 'fontsize', 11, 'fontweight', 'bold');
    ylabel('A_P ({\mu}mol kg^{-1})', 'fontsize', 11, 'fontweight', 'bold');
    grid on;
    legend(h, sprintf("%.1f", glob_b2(1))+"S_A", sprintf("%.1f", glob_b(2))+"S_A + "+sprintf("%.0f", glob_b(1)), 'location', 'northwest');
    xlim([20, 40]);
    ylim([1400, 2800]);
    set(gcf, 'name', mfilename+" 2");

    if SAVE_FIGS
        filename = "sec34_all_fits";
        folder = "../auto_saved/";
        if PER_UNIT_MASS, folder = folder + "per_unit_mass/"; else folder = folder + "total_amount/"; end
        saveas(gcf, folder + filename + ".png");
        print(gcf, folder + filename + ".eps", '-depsc');
    end








    %% Total dissolution rate
    strct = data.accum_data(ALKP, 1);
    strct_alkc_avg = data.accum_data(ALKC_avg, 1);
    strct_alkc_lin = data.accum_data(ALKC_lin, 1);
    strct_alkc_reg = data.accum_data(ALKC_reg, 1);
    smth_fac = 1; % factor for smoothing, number of points on either side for sliding window

    colors = brewermap(9, 'Paired');
    colors = colors([1, 3, 7], :);
    lw = 1.8;

    [ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
    set(gcf, 'name', mfilename+" 3");
    if PER_UNIT_MASS
        xlim(ax1, [-.5 .5])
    else
        xlim(ax0, [-130, 230])
        xlim(ax1, [-.17, .12])
    end
    xlim(ax2, xlim(ax1))
    xticks(ax1, xticks(ax2))

    x = strct.gamma;
    y = [sum(strct.E, 1) + sum(strct.M, 1) - (sum(strct_alkc_avg.E, 1));
        sum(strct.E, 1) + sum(strct.M, 1) - (sum(strct_alkc_lin.E, 1));
        sum(strct.E, 1) + sum(strct.M, 1) - (sum(strct_alkc_reg.E, 1))];
    y = -0.5 * y; % accumulation to dissolution rate
    y_smth = moving_average(y, smth_fac); % smoothen curves

    y_without_bdy = -0.5 * [sum(strct.E(2:4,:), 1) + sum(strct.M, 1) - (sum(strct_alkc_avg.E(2:4,:), 1) + sum(strct_alkc_avg.M, 1));
                    sum(strct.E(2:4,:), 1) + sum(strct.M, 1) - (sum(strct_alkc_lin.E(2:4,:), 1) + sum(strct_alkc_lin.M, 1));
                    sum(strct.E(2:4,:), 1) + sum(strct.M, 1) - (sum(strct_alkc_reg.E(2:4,:), 1) + sum(strct_alkc_reg.M, 1))];

    y_bdy = -0.5 * [strct.E(1,:) - strct_alkc_avg.E(1,:);
                    strct.E(1,:) - strct_alkc_lin.E(1,:);
                    strct.E(1,:) - strct_alkc_reg.E(1,:)];

    y_swr = -0.5 * [data.accum_data(ALKP, 2).E(1,:) - data.accum_data(ALKC_avg, 2).E(1,:);
                    data.accum_data(ALKP, 2).E(1,:) - data.accum_data(ALKC_lin, 2).E(1,:);
                    data.accum_data(ALKP, 2).E(1,:) - data.accum_data(ALKC_reg, 2).E(1,:)];
    
    y_surf = y_bdy - y_swr; % dissolution rate due to surface fluxes (excluding SWR)

    
    % Density plot
    fill_y = [data.sulpis.gamma, fliplr(data.sulpis.gamma)];
    fill_x = [data.sulpis.d_rate + data.sulpis.d_rate_std, fliplr(data.sulpis.d_rate - data.sulpis.d_rate_std)];
    fill(ax0, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
    plot(ax0, data.sulpis.d_rate, data.sulpis.gamma, "Color", [.3 .3 .3], "LineWidth", 1.2);
    for i = 1:size(y_smth, 1)
        plot(ax0, y_smth(i,:), x, "Color", colors(i,:), "Linewidth", lw);
    end
%     plot(ax0, y_without_bdy_smth, x, "Color", [.7 .5 .7], "Linewidth", lw);
    legend(ax0, ["", "", "Alk^*", "WMT (avg)", "WMT (lin)", "WMT (reg)"], "fontsize", 9, "numcolumns", 5, "position", [0.171899484536082,0.880952380952381,0.658928571428571,0.053968253968254]);
    
    % Depth plot
    sulpis_depth_d_rate = data.sulpis.d_rate;
    sulpis_depth_d_rate_std = data.sulpis.d_rate_std;
    if PER_UNIT_MASS
        m = mass_for_densities(x);
        y = y .* m;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
        y_without_bdy = y_without_bdy .* m;
        y_swr = y_swr .* m;
        y_surf = y_surf .* m;
        sulpis_depth_d_rate =  sulpis_depth_d_rate .*  data.sulpis.mass;
        sulpis_depth_d_rate_std =  sulpis_depth_d_rate_std .*  data.sulpis.mass;
    else
        y = y * 0.05;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
        y_without_bdy = y_without_bdy * 0.05;
        y_swr = y_swr * 0.05;
        y_surf = y_surf * 0.05;
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

%     w = exp(-zt/50); % exponentially decaying function with depth
%     y = y_without_bdy + sum(y_swr, 'omitnan') * (w/sum(w));  % divide dissolution from SWR over water column with exponential decay    
    y = y_without_bdy + sum(y_swr, 2, 'omitnan') * O03_param(zt);  % divide dissolution from SWR over water column with exponential decay
    y(:,1) = y(:,1) + sum(y_surf, 2, 'omitnan'); % add surface fluxes to top point
    
    if PER_UNIT_MASS
        m = mass_dz10();
        y = y ./ m; % now divide by mass of each depth layer to return to units of umol/kg/yr
        ys = ys ./ m;
        ys_min = ys_min ./ m;
        ys_max = ys_max ./ m;
        y_without_bdy = y_without_bdy ./ m;
    else
        y = y ./ 10; % divide by dz
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
    
    for i = 1:size(y, 1)
        plot(ax1, y(i,:), zt*1e-3, "Color", colors(i,:), "LineWidth", lw);
        plot(ax2, y(i,:), zt*1e-3, "Color", colors(i,:), "LineWidth", lw);
    end

    if SAVE_FIGS
        filename = "sec34_totals";
        folder = "../auto_saved/";
        if PER_UNIT_MASS, folder = folder + "per_unit_mass/"; else folder = folder + "total_amount/"; end
        saveas(gcf, folder + filename + ".png");
        print(gcf, folder + filename + ".eps", '-depsc');
    end

end

