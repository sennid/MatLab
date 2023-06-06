function pl = convergenceFunc(fn,f,a,b,n, convType)

Video1(1:n) = struct('cdata', [], 'colormap', []);
x = linspace(a, b, 1000);
fig_id = figure;
fig_id.Position = [300 100 900 550];
for i = 1:n

    hold on
    grid on
    plot(x, f(x), 'r-');
    p_id = plot(x, fn(i, x), 'b-');
    
    %m = zeros(1, i);
    %for j = 1:i
    %    m(j) = abs(fn(i, x(j))-f(x(j)));
    %end
    sup = max(abs(fn(i, x)-f(x)));
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 14);
    str = ['n = ', num2str(i), '; sup = ', num2str(sup)];
    ylim([0, 1]);
    title(str)
    set(gca, 'Fontsize', 14)
    hold off
 
    Video1(i) = getframe(gcf);
    delete(p_id)
end

vid0bj = VideoWriter('myVideo.avi');
vid0bj.FrameRate = 5;
vid0bj.Quality = 100;
open(vid0bj);
writeVideo(vid0bj, Video1);
close(vid0bj);

pl = 0;
end