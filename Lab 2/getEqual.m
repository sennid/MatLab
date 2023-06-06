function getEqual(f,g,t0,t1,N)
    h = 0.01;
    t = t0:h:t1;
    tt = [t(1:N-1),t1];
    flag = 1;
    x = f(t);
    y = g(t);
    Video1 = struct('cdata',[],'colormap',[]);
    i = 0;
    mp1 = sqrt(max(diff(f(linspace(t0,t1,N))).^2+diff(g(linspace(t0,t1,N))).^2));
    while flag
        i = i+1;
        px = f(tt);
        py = g(tt);
        x1 = diff(px);
        y1 = diff(py);
        [M, Mpos] = max(x1.^2+y1.^2);
        if (abs(sqrt(M)-sqrt(x1.^2+y1.^2))<0.01)|(Mpos==1)
            flag = 0;
        end
        tt(Mpos) = tt(Mpos)+h;
        hold on;
        grid on;
        p_id = plot(x,y,'k-',px,py,'b*');
        set(gca,'FontSize',14);
        mp = sqrt(M);
        s = ['pdist = ', num2str(mp), '; tdist = ', num2str(mp1)];
        title(s);
        axis equal;
        hold off;
            
        Video1(i) = getframe(gcf);
        if (flag~=0)
            delete(p_id);
        end
    end
    vidObj = VideoWriter('Video6.avi');
    vidObj.FrameRate = 500000;
    vidObj.Quality = 100;
    open(vidObj);
    writeVideo(vidObj,Video1);
    close(vidObj);
end
            