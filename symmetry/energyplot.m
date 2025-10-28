function y = energyplot(energy)

n = length(energy);
figure(1)
semilogy(1:n,energy/sum(energy),'o-','LineWidth',1.5);
hold on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
yticks = get(gca, 'ytick'); % 获取当前y轴刻度值
% newLabels = arrayfun(@(x) sprintf('%.2f%%', x*100), yticks, 'UniformOutput', false); % 将刻度值转换为百分比形式
% set(gca, 'yticklabel', newLabels); % 设置新的y轴刻度标签
ylabel('Energy');
xlabel('k');
set(gcf, 'Units', 'centimeters', 'Position', [15,10,14,10]);
set(gcf,'Color',[1 1 1]);
grid on

figure(2)
semilogy(1:n,cumsum(energy)/sum(energy),'o-','LineWidth',1.5);

ylim([0.4,1]);

set(gca, 'FontName', 'Times New Roman', 'FontSize',12);
ylabel('Total Energy');
xlabel('k');
yticks = get(gca, 'ytick'); % 获取当前y轴刻度值
newLabels = arrayfun(@(x) sprintf('%.2f%%', x*100), yticks, 'UniformOutput', false); % 将刻度值转换为百分比形式
set(gca, 'yticklabel', newLabels); % 设置新的y轴刻度标签
set(gcf, 'Units', 'centimeters', 'Position', [30,10,14,10]);
set(gcf,'Color',[1 1 1]);
ylabel('Total Energy');
xlabel('k');
grid on

end