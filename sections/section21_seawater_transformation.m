function section21_seawater_transformation(data)
load('constants.mat');

% Seawater transformation
x = data.wmt_data(ALK, 1).gamma;
y = data.wmt_data(ALK, 1);
figure; hold on;
plot(x, 0*x, '--k');

x = data.wmt_data(ALK, 1).gamma;
y = data.wmt_data(ALK, 1);
plot(x, 1e-9 * y.G(1,:), LineWidth=2, Color=FLUX_COLORS(1,:))
plot(x, 1e-9 * (y.G(2,:) + y.G(3,:)), LineWidth=2, Color=FLUX_COLORS(2,:))
plot(x, 1e-9 * y.G(4,:), LineWidth=2, Color=FLUX_COLORS(3,:))
grid on;
xlabel('Neutral density (kg m^{-3})', FontSize=12, FontWeight='bold')
ylabel('Transformation (mass Sv)', FontSize=12, FontWeight='bold')
legend('', 'G^{bdy}', 'G^{meso}', 'G^{iso}', FontSize=9, Location='southeast')
xlim([20, 29])
set(gcf, 'name', mfilename);

end