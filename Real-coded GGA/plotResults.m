function plotResults(out, specs, params, elapsedTime)
currentDateTime = datetime('now');
dateTimeStr = string(currentDateTime, 'yyyyMMdd_HHmmSS');
modalityLabel = getModalityLabel(specs);

% Convergence Plot with Subplots
figure('Name', 'GA Convergence Analysis', 'Position', [100, 100, 1200, 500]);

% Subplot 1: Best Cost Convergence
subplot(1, 2, 1);
semilogy(out.bestcost, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]); 
hold on;
grid on;
xlabel('Iterations');
ylabel('Best Cost (logarithmic)');
title(sprintf('Best Cost Convergence - %d Cameras', specs.Cams));
text(0.6*params.MaxIt, max(out.bestcost)*0.5, ...
    sprintf('Final Cost: %.4f\nTime: %.1fsec', out.bestsol.Cost, elapsedTime), ...
    'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;

% Subplot 2: Average Cost Evolution
subplot(1, 2, 2);
semilogy(out.avgcost, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'Population Average'); 
hold on;
semilogy(out.topTenAvgCost, 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'Top 10 Average');
semilogy(out.bestcost, 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410], 'LineStyle', '--', 'DisplayName', 'Best Cost');
grid on;
xlabel('Iterations');
ylabel('Cost (logarithmic)');
title('Population Cost Evolution');
legend('Location', 'best');
hold off;

if specs.warmStart
    warmStr = 'Warm-Started';
else
    warmStr = 'Random Population';
end
sgtitle({sprintf('GA Convergence Analysis - %d Cameras (%s)', specs.Cams, warmStr), modalityLabel}, ...
    'FontSize', 14, 'FontWeight', 'bold');

plotFilename = sprintf('%dCams_Run_%s_convergence.png', specs.Cams, dateTimeStr);
saveas(gcf, plotFilename);
fprintf('Convergence plot saved to: %s\n', plotFilename);

% Camera Visualisation
cameras = cell(specs.Cams, 1);
for i = 1:specs.Cams
    chromStart = (i-1)*6 + 1;
    chromEnd = i*6;
    pos = out.bestsol.Chromosome(chromStart:chromStart+2);
    ori = out.bestsol.Chromosome(chromEnd-2:chromEnd);
    T = se3(eul2rotm(ori, "XYZ"), pos);
    cameras{i} = CentralCamera(name="cam"+i, pose=T);
end

figure('Name', 'Camera Configuration', 'Position', [950, 100, 800, 600]);
hold on;

% Plot cameras
for i = 1:specs.Cams
    cameras{i}.plot_camera('label', scale = 0.5)
end

axis('equal')
grid('on')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title({sprintf('Optimal Camera Placement - Cost: %.4f', out.bestsol.Cost), modalityLabel})
view(45, 30)
hold off

cameraPlotFilename = sprintf('%dCams_Run_%s_cameras.png', specs.Cams, dateTimeStr);
saveas(gcf, cameraPlotFilename); %saves as png
fprintf('Camera plot saved to: %s\n', cameraPlotFilename);