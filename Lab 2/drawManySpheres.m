function F = drawManySpheres(alphas, colors, edges, trans)
    a = -1.01;
    b = 1.01;
    n = 100;
    lev = 1;
    [xMat,yMat,zMat] = meshgrid(linspace(a,b,n),linspace(a,b,n),linspace(a,b,n));
    for i=1:length(alphas)
        if isinf(alphas(i))
            rMat = max(max(abs(xMat), abs(yMat)),abs(zMat));
        else
            rMat = (abs(xMat).^alphas(i) + abs(yMat).^alphas(i) + abs(zMat).^alphas(i)).^(1/alphas(i));
        end
        hold on;
        grid on;
        p = patch(isosurface(xMat,yMat,zMat,rMat,lev));
        p.FaceColor = colors(i);
        p.EdgeColor = edges(i);
        alpha(trans(i));
    end
    xlabel('$x$','Interpreter','latex','FontSize', 14);
    ylabel('$y$','Interpreter','latex','FontSize', 14);
    zlabel('$z$','Interpreter','latex','FontSize', 14);
    view(3)
    camlight
    hold off;
end