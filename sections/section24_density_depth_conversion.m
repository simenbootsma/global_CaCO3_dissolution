function section24_density_depth_conversion(~)
load('constants.mat');
%% Transformation matrix
% [Z, Gamma, A] = compute_matrix();
% save('data\depth_density_matrix_dz10.mat', 'Z', 'Gamma', 'A');
load('depth_density_matrix_dz10.mat', 'Z', 'Gamma', 'A');
zt = sort(unique(Z(:)));
gamma_t = sort(unique(Gamma(:)));

figure;
imagesc(zt*1e-3, gamma_t, A);
grid on 
xlabel('Depth (km)');
ylabel('Neutral density (kg m^{-3})');
caxis([0 .5]);
cb = colorbar;
cb.Label.String = 'Volume fraction (-)';
colors = brewermap(100, '*Spectral');
colors(1,:) = 0;
colormap(colors);
caxis([1e-6 1]);
set(gca, 'colorscale', 'log');
set(gca, 'ydir', 'normal');
pbaspect([2 1 1])
xlim([0 5.5]);
set(gcf, 'name', mfilename + " 1")


%% Distributions
selection = [22, 26.5, 27, 27.5 27.8, 28, 28.2, 28.4 28.9];
figure; hold on;
pbaspect([2 1 1])
colors = brewermap(length(selection), 'YlOrRd');
colors(1:2,:) = colors(1:2,:) * 0.9;
for i = 1:length(selection)
    idx = find(gamma_t >= selection(i), 1);
    y = A(idx, :);
    plot(zt*1e-3, y, '-', 'color', colors(i,:), 'linewidth', 1.5, 'displayname', "\gamma = "+sprintf("%.1f", selection(i))+" kg m^{-3}");
end
grid on;
xlabel('Depth (km)', 'fontweight', 'bold');
ylabel('Volume fraction (-)','fontweight', 'bold');
legend('Numcolumns', 3);
set(gcf, 'name', mfilename + " 2")

end




function [Z, Gamma, A] = compute_matrix()
    load('WOA18_NeutralDensity3D.mat', 'gamma');
    load('WOA18_Volume.mat', 'vol');
    load('WOA18_Coordinates.mat', 'z');
    vol = repmat(vol, [1,1,1,12]);
    vol2d = vol(:,:,1);
    
    equal_z = 10:10:5500;
    equal_gamma = nan(360,180,length(equal_z));
    
    for i = 1:size(gamma,1)
        for j = 1:size(gamma,2)
            gamma_column = reshape(gamma(i,j,:), [1,length(z)]);
            intp_gamma = interp1(z, gamma_column, equal_z);
            equal_gamma(i,j,:) = intp_gamma;
        end
    end
    clear vol gamma;
    
    gamma_edges = 14:0.05:30;
    gamma_t = gamma_edges(1:end-1) + diff(gamma_edges)/2;
    
    A = zeros(length(gamma_t), length(equal_z));
    [Z, Gamma] = meshgrid(equal_z, gamma_t);
    for i = 1:size(A,2)
        tic;
        egi = equal_gamma(:,:,i);
        for j = 1:size(A,1)
            cnts = sum((egi > gamma_edges(j)) & (egi <= gamma_edges(j+1)), 3);
            A(j,i) = sum(vol2d.*cnts, 'all');
        end
        fprintf("%d/%d: %.1f seconds\n", [i, size(A,2), toc]);
    end
    A = A ./ sum(A, 2, 'omitnan');
    A(isnan(A)) = 0;
end