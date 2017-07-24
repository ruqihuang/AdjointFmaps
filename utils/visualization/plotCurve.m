function h = plotCurve(x, y)

h = plot(x, y);
set(h,'linewidth',3);
grid on;


ylim([0 100]);
ylabel('% of Correspondences', 'FontSize',25);
xlabel('% of Geodesic Error', 'FontSize',25);
% leg = legend('Ground Truth Fmap', 'Fmap', 'Location','southeast');
% set(leg,'FontSize',25);
hold off;
end