function [ax_left, ax_right1, ax_right2] = setup_profiles_plot(PER_UNIT_MASS)

%% Density plot
figure;
ax_left = subtightplot(2, 2, [1 3], 0.001, [0.128, 0.15], 0.1);
hold on;
plot(ax_left, [0, 0], [20, 29], '--k');
xlabel(ax_left, "r$_{CaCO_3}$ (Tmol $\Delta\gamma^{-1}$ yr$^{-1}$)", FontSize=11, FontWeight='bold', Interpreter='latex');
ylabel(ax_left, "Neutral density (kg m^{-3})", FontSize=11, FontWeight='bold');
grid(ax_left);
set(gca, 'ydir', 'reverse')
ylim(ax_left, [20 29])
pbaspect(ax_left, [1 2 1]);


%% Depth plot
ax_right1 = subtightplot(2, 2, 2, 0.01, [0.128, 0.15], 0.1);
hold on;
plot(ax_right1, [0, 0], [0, 1], '--k');

ylabel(ax_right1, "Depth (km)", FontSize=11, FontWeight='bold', Units='normalized', Position=[-.18,0,0]);
grid(ax_right1);
set(gca, 'ydir', 'reverse')
ylim(ax_right1, [0 1])
xlim(ax_right1, xlim(ax_left));
xticklabels(ax_right1, []);

ytickformat(ax_right1, '%.1f');
pbaspect(ax_right1, [1 1 1]);

ax_right2 = subtightplot(2, 2, 4, 0.01, [0.128, 0.15], 0.1);
hold on;
plot(ax_right2, [0, 0], [1, 6], '--k');

xlabel(ax_right2, "r$_{CaCO_3}$ (Tmol m$^{-1}$ yr$^{-1}$)", FontSize=11, FontWeight='bold', Interpreter='latex');
grid(ax_right2);
set(ax_right2, 'ydir', 'reverse')
ylim(ax_right2, [1.01 6])
pbaspect(ax_right2, [1 1 1]);
xlim(ax_right2, xlim(ax_left));
xticks(ax_right1, xticks(ax_right2));
ytickformat(ax_right2, '%.1f');

if PER_UNIT_MASS
    xlabel(ax_left, "r$_{CaCO_3}$ ($\mu$mol kg$^{-1}$ yr$^{-1}$)", FontSize=11, FontWeight='bold', Interpreter='latex');
    xlabel(ax_right2, "r$_{CaCO_3}$ ($\mu$mol kg$^{-1}$ yr$^{-1}$)", FontSize=11, FontWeight='bold', Interpreter='latex');
end

end