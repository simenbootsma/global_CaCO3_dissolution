function section32_dissolution_profiles_wmt(data)
    load('constants.mat');
    strct = data.accum_data(ALKP, 1);
    lw = 1.5;
    smth_fac = 3; % factor for smoothing, number of points on either side for sliding window
    
    %% 4-panel with WMT curves

    % WMT of Alk
    figure; hold on;
    x = data.wmt_data(ALK, 1).gamma;
    y = data.wmt_data(ALK, 1);
    y.E = moving_average(y.E, smth_fac);
    y.M = moving_average(y.M, smth_fac);

    plot(x, 0*x, '--k');
    plot(x, y.E(1,:), LineWidth=1.5, Color=FLUX_COLORS(1,:))
    plot(x, y.E(2,:) + y.E(3,:), LineWidth=1.5, Color=FLUX_COLORS(2,:))
    plot(x, y.E(4,:), LineWidth=1.5, Color=FLUX_COLORS(3,:))
    plot(x, y.M(1,:) + y.M(2,:), LineWidth=1.5, Color=FLUX_COLORS(4,:))
    plot(x, y.M(3,:), LineWidth=1.5, Color=FLUX_COLORS(5,:))
    grid on;
    xlabel('Neutral density (kg m^{-3})', FontSize=12, FontWeight='bold')
    ylabel('Diapycnal transport of A_T (Tmol yr^{-1})', FontSize=12, FontWeight='bold')
    legend('', 'E^{bdy}', 'E^{meso}', 'E^{iso}', 'M^{meso}', 'M^{iso}', FontSize=9, NumColumns=5, Position=[0.18,0.935,0.7,0.05])
    xlim([20, 29])
    set(gcf, 'name', mfilename+" 1");



    % WMT of NO3
    x = data.wmt_data(NIT, 1).gamma;
    y = data.wmt_data(NIT, 1);
    y.E = moving_average(y.E, smth_fac);
    y.M = moving_average(y.M, smth_fac);

    figure; hold on;
    plot(x, 0*x, '--k');
    plot(x, y.E(1,:), LineWidth=1.5, Color=FLUX_COLORS(1,:))
    plot(x, (y.E(2,:) + y.E(3,:)), LineWidth=1.5, Color=FLUX_COLORS(2,:))
    plot(x, y.E(4,:), LineWidth=1.5, Color=FLUX_COLORS(3,:))
    plot(x, (y.M(1,:) + y.M(2,:)), LineWidth=1.5, Color=FLUX_COLORS(4,:))
    plot(x, y.M(3,:), LineWidth=1.5, Color=FLUX_COLORS(5,:))
    grid on;
    xlabel('Neutral density (kg m^{-3})', FontSize=12, FontWeight='bold')
    ylabel('Diapycnal transport of NO_3 (Tmol yr^{-1})', FontSize=12, FontWeight='bold')
    legend('', 'E^{bdy}', 'E^{meso}', 'E^{iso}', 'M^{meso}', 'M^{iso}', FontSize=9, NumColumns=5, Position=[0.18,0.935,0.7,0.05])
    xlim([20, 29])
    ylim([-25, 15])
    set(gcf, 'name', mfilename+" 2");


    % WMT of AlkP
    x = data.wmt_data(ALKP, 1).gamma;
    y = data.wmt_data(ALKP, 1);
    y.E = moving_average(y.E, smth_fac);
    y.M = moving_average(y.M, smth_fac);

    figure; hold on;
    plot(x, 0*x, '--k');
    plot(x, y.E(1,:), LineWidth=1.5, Color=FLUX_COLORS(1,:))
    plot(x, (y.E(2,:) + y.E(3,:)), LineWidth=1.5, Color=FLUX_COLORS(2,:))
    plot(x, y.E(4,:), LineWidth=1.5, Color=FLUX_COLORS(3,:))
    plot(x, (y.M(1,:) + y.M(2,:)), LineWidth=1.5, Color=FLUX_COLORS(4,:))
    plot(x, y.M(3,:), LineWidth=1.5, Color=FLUX_COLORS(5,:))
    grid on;
    xlabel('Neutral density (kg m^{-3})', FontSize=12, FontWeight='bold')
    ylabel('Diapycnal transport of A_P (Tmol yr^{-1})', FontSize=12, FontWeight='bold')
    legend('', 'E^{bdy}', 'E^{meso}', 'E^{iso}', 'M^{meso}', 'M^{iso}', FontSize=9, NumColumns=5, Position=[0.18,0.935,0.7,0.05])
    xlim([20, 29])
    set(gcf, 'name', mfilename+" 3");



    % WMF of AlkP
    smth_fac2 = 7;
    x = data.accum_data(ALKP, 1).gamma;
    y = data.accum_data(ALKP, 1);
    y.E = moving_average(y.E, smth_fac2);
    y.M = moving_average(y.M, smth_fac2);

    figure; hold on;
    plot(x, 0*x, '--k');
    plot(x, y.E(1,:), LineWidth=1.5, Color=FLUX_COLORS(1,:))
    plot(x, (y.E(2,:) + y.E(3,:)), LineWidth=1.5, Color=FLUX_COLORS(2,:))
    plot(x, y.E(4,:), LineWidth=1.5, Color=FLUX_COLORS(3,:))
    plot(x, (y.M(1,:) + y.M(2,:)), LineWidth=1.5, Color=FLUX_COLORS(4,:))
    plot(x, y.M(3,:), LineWidth=1.5, Color=FLUX_COLORS(5,:))
    grid on;
    xlabel('Neutral density (kg m^{-3})', FontSize=12, FontWeight='bold')
    if PER_UNIT_MASS
        ylabel('Accumulation of A_P ({\mu}mol kg^{-1} yr^{-1})', FontSize=12, FontWeight='bold')
    else
        ylabel('Accumulation of A_P (Tmol \Delta\gamma^{-1} yr^{-1})', FontSize=12, FontWeight='bold')
    end
    legend('', '{\Delta}E^{bdy}', '{\Delta}E^{meso}', '{\Delta}E^{iso}', '{\Delta}M^{meso}', '{\Delta}M^{iso}', FontSize=9, NumColumns=5, Position=[0.13,0.943,0.78,0.05])
    xlim([20, 29])
    set(gcf, 'name', mfilename+" 4");








    %% All fluxes 
    [ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
    if PER_UNIT_MASS
        xlim(ax0, [-1000 1000]);
        xlim(ax1, [-30 50]);
    else
        xlim(ax0, [-7000, 7000])
        xlim(ax1, [-12 12])
    end
    xlim(ax2, xlim(ax1));
    xticks(ax1, xticks(ax2));
    set(gcf, 'name', mfilename+" 5");

    x = strct.gamma;
    y = [strct.E(1,:);
        strct.E(2,:) + strct.E(3,:);
        strct.E(4,:);
        strct.M(1,:)+strct.M(2,:);
        strct.M(3,:)];
    y = -0.5 * y; % accumulation to dissolution rate
    y_swr = -0.5 * data.accum_data(ALKP, 2).E(1,:);  % dissolution rate due to short-wave radiation
    y_surf = y(1,:) - y_swr; % dissolution rate due to surface fluxes (excluding SWR)
    y_smth = moving_average(y, smth_fac); % smoothen curves
    
    % Density plot
    for i = 1:size(y,1)
        plot(ax0, y_smth(i,:), x, "Color", FLUX_COLORS(i,:), "Linewidth", lw);
    end
    legend(ax0, ["", "E^{bdy}", "E^{meso}", "E^{iso}", "M^{meso}", "M^{iso}"], "fontsize", 9, "NumColumns", 5, "Position", [0.155232817869416,0.877777777777777,0.697023809523809,0.053968253968254]);
    
    % Depth plot
    if PER_UNIT_MASS
        m = mass_for_densities(x);
        y = y .* m;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
        y_swr = y_swr .* m;
        y_surf = y_surf .* m;
    else
        y = y * 0.05;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
        y_swr = y_swr * 0.05;
        y_surf = y_surf * 0.05;
    end
    load('depth_density_matrix_dz10.mat', 'Z', 'Gamma', 'A');
    x(isnan(x)) = 0;
    y(isnan(y)) = 0;
    zt = sort(unique(Z(:)))';
    B = interp2(Z, Gamma, A, zt, x', "linear", 0);
    y = y * B;

    y(1,:) = sum(y_swr, 'omitnan') * O03_param(zt);
    y(1,1) = y(1,1) + sum(y_surf, 'omitnan'); % add surface fluxes to top point
    
    if PER_UNIT_MASS
        m = mass_dz10();
        y = y ./ m; % now divide by mass of each depth layer to return to units of umol/kg/yr
    else
        y = y / 10; % divide by dz
    end

    for i = 1:size(y,1)
        plot(ax1, y(i,:), zt*1e-3, "Color", FLUX_COLORS(i,:), "LineWidth", lw);
        plot(ax2, y(i,:), zt*1e-3, "Color", FLUX_COLORS(i,:), "LineWidth", lw);
    end




    %% Total dissolution rate
    colors = brewermap(2, '*Pastel1') * .9;
    lw = 1.8;

    [ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
    if PER_UNIT_MASS
        xlim(ax0, [-1000 1000]);
        xlim(ax1, [-30 50]);
    else
        xlim(ax0, [-7000, 7000])
        xlim(ax1, [-12 12])
    end
    xlim(ax2, xlim(ax1));
    xticks(ax1, xticks(ax2));

    set(gcf, 'name', mfilename+" 6");

    x = strct.gamma;
    y = sum(strct.E, 1) + sum(strct.M, 1);
    y = -0.5 * y; % accumulation to dissolution rate
    y_smth = moving_average(y, smth_fac); % smoothen curves
    y_without_bdy = -0.5 * (sum(strct.E(2:4, :), 1) + sum(strct.M, 1));
    y_without_bdy_smth = moving_average(y_without_bdy, smth_fac);
    y_bdy = -0.5 * data.accum_data(ALKP, 1).E(1,:);  % dissolution rate due to boundary fluxes
    y_swr = -0.5 * data.accum_data(ALKP, 2).E(1,:);  % dissolution rate due to short-wave radiation
    y_surf = y_bdy - y_swr; % dissolution rate due to surface fluxes (excluding SWR)
    
    % Density plot
    fill_y = [data.sulpis.gamma, fliplr(data.sulpis.gamma)];
    fill_x = [data.sulpis.d_rate + data.sulpis.d_rate_std, fliplr(data.sulpis.d_rate - data.sulpis.d_rate_std)];
    fill(ax0, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
    plot(ax0, data.sulpis.d_rate, data.sulpis.gamma, "Color", [.3 .3 .3], "LineWidth", 1.2);

    plot(ax0, y_smth, x, "Color", colors(1,:), "Linewidth", lw);
    plot(ax0, y_without_bdy_smth, x, "Color", colors(2,:), "Linewidth", lw);
    legend(ax0, ["", "", "Alk^*", "WMT (with E^{bdy})", "WMT (without E^{bdy})"], "fontsize", 9, "numcolumns", 5, "position", [0.171899484536082,0.880952380952381,0.658928571428571,0.053968253968254]);
    
    % Depth plot
    if PER_UNIT_MASS
        m = mass_for_densities(x);
        y = y .* m;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
        y_without_bdy = y_without_bdy .* m;
        y_swr = y_swr .* m;
        y_surf = y_surf .* m;
        data.sulpis.d_rate =  data.sulpis.d_rate .*  data.sulpis.mass;
        data.sulpis.d_rate_std =  data.sulpis.d_rate_std .*  data.sulpis.mass;
    else
        y = y * 0.05;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
        y_without_bdy = y_without_bdy * 0.05;
        y_swr = y_swr * 0.05;
        y_surf = y_surf * 0.05;
        data.sulpis.d_rate =  data.sulpis.d_rate * 0.05;
        data.sulpis.d_rate_std =  data.sulpis.d_rate_std * 0.05;
    end
    load('depth_density_matrix_dz10.mat', 'Z', 'Gamma', 'A');
    data.sulpis.gamma(isnan(data.sulpis.gamma)) = 0;
    data.sulpis.d_rate(isnan(data.sulpis.d_rate)) = 0;
    zt = sort(unique(Z(:)))';
    B_sulpis = interp2(Z, Gamma, A, zt, data.sulpis.gamma', "linear", 0);
    ys = data.sulpis.d_rate * B_sulpis;
    ys_min = (data.sulpis.d_rate - data.sulpis.d_rate_std) * B_sulpis;
    ys_max = (data.sulpis.d_rate + data.sulpis.d_rate_std) * B_sulpis;

    x(isnan(x)) = 0;
    y(isnan(y)) = 0;
    y_without_bdy(isnan(y_without_bdy)) = 0;
    B = interp2(Z, Gamma, A, zt, x', "linear", 0);
    y_without_bdy = y_without_bdy * B;

    y = y_without_bdy + sum(y_swr, 'omitnan') * O03_param(zt);  % divide dissolution from SWR over water column with exponential decay
    y(1) = y(1) + sum(y_surf, 'omitnan'); % add surface fluxes to top point
    
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

    plot(ax1, y, zt*1e-3, "Color", colors(1,:), "LineWidth", lw);
    plot(ax1, y_without_bdy, zt*1e-3, "Color", colors(2,:), "LineWidth", lw);
    plot(ax2, y, zt*1e-3, "Color", colors(1,:), "LineWidth", lw);
    plot(ax2, y_without_bdy, zt*1e-3, "Color", colors(2,:), "LineWidth", lw);
    
end

