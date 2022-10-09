function section33_error_compensation(data)
    load('constants.mat');
    strct = data.accum_data(ALKP, 1);
    strct_sal = data.accum_data(SAL, 1);
    lw = 1.5;
    smth_fac = 3; % factor for smoothing, number of points on either side for sliding window

    %% WMF of Sal
    smth_fac2 = 7;
    x = data.accum_data(SAL, 1).gamma;
    y = data.accum_data(SAL, 1);
    y.E = moving_average(y.E, smth_fac2);
    y.M = moving_average(y.M, smth_fac2);

    figure; hold on;
    plot(x, 0*x, '--k');
    plot(x, y.E(1,:), LineWidth=1.5, Color=FLUX_COLORS(1,:))
    plot(x, (y.M(1,:) + y.M(2,:)), LineWidth=1.5, Color=FLUX_COLORS(4,:))
    plot(x, (y.E(2,:) + y.E(3,:)), LineWidth=1.5, Color=FLUX_COLORS(2,:))
    plot(x, y.M(3,:), LineWidth=1.5, Color=FLUX_COLORS(5,:))
    plot(x, y.E(4,:), LineWidth=1.5, Color=FLUX_COLORS(3,:))
    plot(x, sum(y.E, 1)+sum(y.M, 1), '-', LineWidth=2, Color=[.3 .3 .3])
    grid on;
    xlabel('Neutral density (kg m^{-3})', FontSize=12, FontWeight='bold')
    if PER_UNIT_MASS
        ylabel('Accumulation of salt (g kg^{-1} yr^{-1})', FontSize=12, FontWeight='bold')
    else
        ylabel('Accumulation of salt (Tg \Delta\gamma^{-1} yr^{-1})', FontSize=12, FontWeight='bold')
    end
    legend('', '\DeltaE^{bdy}', '\DeltaM^{meso}', '\DeltaE^{meso}', '\DeltaM^{iso}', '\DeltaE^{iso}', '\DeltaE+\DeltaM', FontSize=9, NumColumns=3, Position=[0.308571428571429,0.793307855151948,0.58,0.099543019854833])

    xlim([20, 29])
    set(gcf, 'name', mfilename+" 4");




    %% All fluxes 
    [ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
    set(gcf, 'name', mfilename+" 1");
    if PER_UNIT_MASS
        xlim(ax1, [-.5 .5])
    else
        xlim(ax0, [-100, 200])
        xlim(ax1, [-.2, .2])
    end
    xlim(ax2, xlim(ax1))
    xticks(ax1, xticks(ax2))

    x = strct.gamma;
    y = [strct.E(1,:) - 66.1 * strct_sal.E(1,:);
        strct.E(2,:) + strct.E(3,:) - 66.1 * (strct_sal.E(2,:) + strct_sal.E(3,:));
        strct.E(4,:) - 66.1 * strct_sal.E(4,:);
        strct.M(1,:)+strct.M(2,:) - 66.1 * (strct_sal.M(1,:) + strct_sal.M(2,:));
        strct.M(3,:) - 66.1 * strct_sal.M(3,:)];

    y = -0.5 * y; % accumulation to dissolution rate
    y_swr = -0.5 * (data.accum_data(ALKP, 2).E(1,:) - 66.1 * data.accum_data(SAL, 2).E(1,:));  % dissolution rate due to short-wave radiation
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
        y = y .* 0.05;  % convert to total dissolution (umol/yr) before transforming into depth coordinates
        y_swr = y_swr .* 0.05;
        y_surf = y_surf .* 0.05;
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
        y = y ./ 10; % divide by dz
    end

    for i = 1:size(y,1)
        plot(ax1, y(i,:), zt*1e-3, "Color", FLUX_COLORS(i,:), "LineWidth", lw);
        plot(ax2, y(i,:), zt*1e-3, "Color", FLUX_COLORS(i,:), "LineWidth", lw);
    end




    %% Total dissolution rate
    colors = brewermap(2, '*Pastel1') * .9;
    lw = 1.8;

    [ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
    set(gcf, 'name', mfilename+" 2");
    if PER_UNIT_MASS
        xlim(ax1, [-.5 .5])
    else
        xlim(ax0, [-100, 200])
        xlim(ax1, [-.17, .12])
    end
    xlim(ax2, xlim(ax1))
    xticks(ax1, xticks(ax2))

    x = strct.gamma;
    y = sum(strct.E, 1) + sum(strct.M, 1) - 66.1 * (sum(strct_sal.E, 1) + sum(strct_sal.M, 1));
    y = -0.5 * y; % accumulation to dissolution rate
    y_smth = moving_average(y, smth_fac); % smoothen curves
    y_without_bdy = -0.5 * (sum(strct.E(2:4, :), 1) + sum(strct.M, 1) - 66.1 * (sum(strct_sal.E(2:4, :), 1) + sum(strct_sal.M, 1)));
    y_without_bdy_smth = moving_average(y_without_bdy, smth_fac);
    y_bdy = -0.5 * (data.accum_data(ALKP, 1).E(1,:) - 66.1 * data.accum_data(SAL, 1).E(1,:));  % dissolution rate due to boundary fluxes
    y_swr = -0.5 * (data.accum_data(ALKP, 2).E(1,:) - 66.1 * data.accum_data(SAL, 2).E(1,:));  % dissolution rate due to short-wave radiation
    y_surf = y_bdy - y_swr; % dissolution rate due to surface fluxes (excluding SWR)
    
    % Density plot
    fill_y = [data.sulpis.gamma, fliplr(data.sulpis.gamma)];
    fill_x = [data.sulpis.d_rate + data.sulpis.d_rate_std, fliplr(data.sulpis.d_rate - data.sulpis.d_rate_std)];
    fill(ax0, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
    plot(ax0, data.sulpis.d_rate, data.sulpis.gamma, "Color", [.3 .3 .3], "LineWidth", 1.2);

%     disp("dissolution (Sulpis): "+string(sum(data.sulpis.d_rate, 'omitnan')*0.05)+" Tmol/yr")

    plot(ax0, y_smth, x, "Color", colors(1,:), "Linewidth", lw);
    plot(ax0, y_without_bdy_smth, x, "Color", colors(2,:), "Linewidth", lw);
    legend(ax0, ["", "", "Alk^*", "WMT (with E^{bdy})", "WMT (without E^{bdy})"], "fontsize", 9, "numcolumns", 5, "position", [0.171899484536082,0.880952380952381,0.658928571428571,0.053968253968254]);
   
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






    %% Total dissolution rate (zoomed in)
    lw = 1.8;

    [ax0, ax1, ax2] = setup_profiles_plot(PER_UNIT_MASS);
    set(gcf, 'name', mfilename+" 3");
    if PER_UNIT_MASS
        xlim(ax1, [-.5 .5])
    else
        xlim(ax0, [-60, 80])
        xlim(ax1, [-.15, .0501])
    end
    xlim(ax2, xlim(ax1))
    xticks(ax1, xticks(ax2))

    x = strct.gamma;
    y = sum(strct.E, 1) + sum(strct.M, 1) - 66.1 * (sum(strct_sal.E, 1) + sum(strct_sal.M, 1));
    y = -0.5 * y; % accumulation to dissolution rate
    y_without_bdy = -0.5 * (sum(strct.E(2:4, :), 1) + sum(strct.M, 1) - 66.1 * (sum(strct_sal.E(2:4, :), 1) + sum(strct_sal.M, 1)));
    y_without_bdy_smth = moving_average(y_without_bdy, smth_fac);
    y_bdy = -0.5 * (data.accum_data(ALKP, 1).E(1,:) - 66.1 * data.accum_data(SAL, 1).E(1,:));  % dissolution rate due to boundary fluxes
    y_swr = -0.5 * (data.accum_data(ALKP, 2).E(1,:) - 66.1 * data.accum_data(SAL, 2).E(1,:));  % dissolution rate due to short-wave radiation
    y_surf = y_bdy - y_swr; % dissolution rate due to surface fluxes (excluding SWR)
    
    % Density plot
%     fill_y = [data.sulpis.gamma, fliplr(data.sulpis.gamma)];
%     fill_x = [data.sulpis.d_rate + data.sulpis.d_rate_std, fliplr(data.sulpis.d_rate - data.sulpis.d_rate_std)];
%     fill(ax0, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
%     plot(ax0, data.sulpis.d_rate, data.sulpis.gamma, "Color", [.3 .3 .3], "LineWidth", 1.2);

    I = 1:length(x);
    yfill = y_without_bdy_smth(I);

%     disp("dissolution: "+string(sum(y(y > 0))*.05)+" Tmol/yr")
%     disp("calcification: "+string(sum(y(y < 0))*.05)+" Tmol/yr")
%     disp("dissolution (no bdy): "+string(sum(yfill(yfill > 0))*.05)+" Tmol/yr")
%     disp("calcification (no bdy): "+string(sum(yfill(yfill < 0))*.05)+" Tmol/yr")

    yfill(yfill > 0) = 0;
    fill(ax0, yfill, x(I), [0 0 1], 'FaceAlpha', .5);
    yfill = y_without_bdy_smth(I);
    yfill(yfill < 0) = 0;
    fill(ax0, yfill, x(I), [1 0 0], 'FaceAlpha', .5);

    plot(ax0, y_without_bdy_smth, x, "Color", colors(2,:), "Linewidth", lw);
%     legend(ax0, ["", "", "Alk^*", "", "", "WMT (without E^{bdy})"], "fontsize", 9, "numcolumns", 5, "position", [0.171899484536082,0.880952380952381,0.658928571428571,0.053968253968254]);

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

%     fill_y = [zt, fliplr(zt)] * 1e-3;
%     fill_x = [ys_max, fliplr(ys_min)];
%     fill(ax1, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
%     plot(ax1, ys, zt*1e-3, "Color", [.3 .3 .3], "LineWidth", 1.2);
%     fill(ax2, fill_x, fill_y, [.4 .4 .4], 'FaceAlpha', .1);
%     plot(ax2, ys, zt*1e-3, "Color", [.3 .3 .3], "LineWidth", 1.2);

%     disp("dissolution: "+string(sum(y(y > 0))*10)+" Tmol/yr")
%     disp("calcification: "+string(sum(y(y < 0))*10)+" Tmol/yr")
%     disp("dissolution (no bdy): "+string(sum(y_without_bdy(y_without_bdy > 0))*10)+" Tmol/yr")
%     disp("calcification (no bdy): "+string(sum(y_without_bdy(y_without_bdy < 0))*10)+" Tmol/yr")

    disp("net dissolution:")
    disp("  A: "+string(10*sum(y(zt<350)))+" Tmol/yr")
    disp("  B: "+string(10*sum(y(zt>=350 & zt<=1020)))+" Tmol/yr")
    disp("  C: "+string(10*sum(y(zt>1020&zt<1530)))+" Tmol/yr")
    disp("  D: "+string(10*sum(y(zt>=1530&zt<=3310)))+" Tmol/yr")
    disp("  E: "+string(10*sum(y(zt>3310)))+" Tmol/yr")

    yfill = [0, y];
    yfill(yfill > 0) = 0;
    fill(ax1, yfill, [zt(1), zt]*1e-3, [0 0 1], 'FaceAlpha', .5);
    fill(ax2, yfill, [zt(1), zt]*1e-3, [0 0 1], 'FaceAlpha', .5);
    yfill = [0, y];
    yfill(yfill < 0) = 0;
    fill(ax1, yfill, [zt(1), zt]*1e-3, [1 0 0], 'FaceAlpha', .5);
    fill(ax2, yfill, [zt(1), zt]*1e-3, [1 0 0], 'FaceAlpha', .5);

    plot(ax1, y, zt*1e-3, "Color", colors(2,:), "LineWidth", lw);
    plot(ax2, y, zt*1e-3, "Color", colors(2,:), "LineWidth", lw);

end

