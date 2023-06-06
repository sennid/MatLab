function flag = makeaplot(x_vec, y_vec)


fig_id = figure;
fig_id.Position = [300 100 1000 650];
hold on
grid on
title('$T(n)$', 'Interpreter', 'latex');
plot(x_vec, y_vec, 'r-', 'Linewidth', 1);
xlabel('$n$', 'Interpreter','latex');
ylabel('$T(n)$', 'Interpreter','latex');
xlim([0, max(x_vec) + 1]);
ylim([min(y_vec) - 1, max(y_vec) + 1]);

hold off

legend({'$T(n)$'}, "Interpreter", "latex", 'FontSize', 15, 'Location', 'northeast');

flag = 0;