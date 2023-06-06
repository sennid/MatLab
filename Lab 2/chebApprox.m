function F = chebApprox(f,n)
    q = 10000;
    TMat = getFunc(n); %матрица коэффициентов полиномов
    web = linspace(-1,1,q);
    y = f(web);

    Video1(1:n+1) = struct('cdata',[],'colormap',[]);
    Sn = zeros(1,q);
    normVec = sqrt(1-0.9999*web.^2);

    for k = 0:n
        Cheb = zeros(1, q);
        for i = 1:n+1
            Cheb = Cheb+TMat(k+1, i)*web.^(i-1);
        end
        Sn = Sn + Cheb.*trapz(web, y.*Cheb./normVec)./trapz(web, Cheb.^2./normVec);

        hold on;
        grid on;


        p_id = plot(web,y,'b-',web,Sn,'g-');
        xlabel('$x$','Interpreter','latex','FontSize', 14);
        ylabel('$f(x)$','Interpreter','latex','FontSize', 14);
        legend({'$f(x)$','$S_n(x)$'},'Interpreter','latex', 'FontSize', 15, 'Location', 'northeast');
        str = num2str(k);
        title(str);
    
        set(gca,'FontSize',14);
        hold off;
            
        Video1(k+1) = getframe(gcf);
        if k~=n
            delete(p_id);
        end
    end

    vidObj = VideoWriter('Video4.avi');
    vidObj.FrameRate = 10;
    vidObj.Quality = 100;
    open(vidObj);
    writeVideo(vidObj,Video1);
    close(vidObj);
end