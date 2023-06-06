%% task12
clc; clear; format compact; close all
f = @(x) exp(-x.^2);

n = 1000;
xl = -100;
xr = 100;

x1Vec = linspace(xl, xr, n);

tic
sT = cumtrapz(x1Vec, f(x1Vec));
tT = toc;

tic  
sR = rectangless(x1Vec, f);
tR = toc;

tic 
sS = simpsonn(x1Vec, f);
tS = toc;

fig_id = figure;
fig_id.Position = [300 100 900 550];
hold on
grid on
title('Integral')
plot(x1Vec,sT,'r-')
plot(x1Vec,sS,'b-')
plot(x1Vec,sR,'g-')
plot([xl xr], [(pi/2)^0.5, (pi/2)^0.5], 'k-')
legend('trapz', 'simpson', 'rectangles')
xlabel('$x$','Interpreter','Latex')
ylabel('$integral$','Interpreter','Latex')
hold off

hVec = 0.001:0.01:1;
m = length(hVec);
rTVec = ones(1,m);
rSVec = ones(1,m);
rRVec = ones(1,m);

for i = 1:m
    x1Vec = xl:hVec(i):xr;
    x2Vec = xl:(hVec(i)/2):xr;
    cumtrapz1 = cumtrapz(x1Vec, f(x1Vec));
    cumtrapz2 = cumtrapz(x2Vec, f(x2Vec));
    simpson1 = simpsonn(x1Vec, f);
    simpson2 = simpsonn(x2Vec, f);
    rec1 = rectangless(x1Vec, f);
    rec2 = rectangless(x2Vec, f);
    rT = (cumtrapz1(end)-cumtrapz2(end));
    rS = (simpson1(end)-simpson2(end)); 
    rR = (rec1(end)-rec2(end)) ;
    rTVec(i) = rT(end);
    rSVec(i) = rS(end);
    rRVec(i) = rR(end);
end

fig_id1 = figure;
fig_id1.Position = [300 100 900 550];
hold on
grid on
title('error')
plot(hVec,rTVec,'r-')
plot(hVec,rSVec,'b-')
plot(hVec,rRVec,'g-')
legend('trapz', 'simpson', 'rectangles')
xlabel('$h$','Interpreter','Latex')
ylabel('$error$','Interpreter','Latex')
hold off