function pl = compareInterp(x,xx,f)

yy1 = interp1(x, f(x),  xx, 'nearest');
yy2 = interp1(x, f(x),  xx, 'linear');
yy3 = interp1(x, f(x),  xx, 'spline');
yy4 = interp1(x, f(x),  xx, 'cubic');

fig_id = figure;
fig_id.Position = [300 100 900 550];
hold on
grid on
title('Interpolation')
plot(xx, f(xx),'o', xx, yy1, '-', xx, yy2, '-', xx, yy3, '-', xx, yy4, '-')
legend('func','nearest', 'linear', 'spline', 'cubic')
xlabel('$x$','Interpreter','Latex')
ylabel('$integral$','Interpreter','Latex')
hold off

legend({'func','nearest', 'linear', 'spline', 'cubic'}, "Interpreter", "latex", 'FontSize', 15, 'Location', 'northeast');

pl = 0;

end