function F = drawSphere(alpha, params)
    n = params.num;
    a = params.a;
    b = params.b;
    lev = params.lev;
    [xMat, yMat, zMat] = meshgrid(linspace(a,b,n),linspace(a,b,n),linspace(a,b,n));
    if isinf(alpha)
        rMat = max(max(abs(xMat), abs(yMat)),abs(zMat));
    else
        rMat = (abs(xMat).^alpha + abs(yMat).^alpha + abs(zMat).^alpha).^(1/alpha);
    end
    figure;
    hold on;
    grid on;
    p = patch(isosurface(xMat,yMat,zMat,rMat,lev));
    p.FaceColor = params.color;
    p.EdgeColor = 'None';
    xlabel('$x$','Interpreter','latex','FontSize', 14);
    ylabel('$y$','Interpreter','latex','FontSize', 14);
    zlabel('$z$','Interpreter','latex','FontSize', 14);
    view(3)
    camlight
    axis equal
    hold off;
end