%% Block 1 (task 1).
clc; clear; format compact; close all

%-----Parameters-------
N = 1000;
step = 1;
S = log(2);
%----End of parameters-----

SVec = ((-ones(1,N)).^(0:N-1))./(1:N);
xVec = 1:N;
PsiVec = 1./(1+xVec);

fig_id = figure;
fig_id.Position = [300 100 1000 650];

hold on 
grid on
plot(xVec, cumsum(SVec)-S, '-r','Linewidth', 1);
plot(xVec, PsiVec, '-b', 'Linewidth', 1);
xlabel('n', 'Interpreter','latex');
ylabel('Sn', 'Interpreter','latex');
hold off

%% Block 2 (task 2).
clc; clear; format compact; close all

%-----Parameters-------
f = @(x) x.^3+x-1-sin(x);
n = 1000;
a = -3;
b = 3;
%----End of parameters-----


xVec = linspace(a,b,n);
fig_id = figure;
fig_id.Position = [300 100 1000 650];

hold on
grid on
plot(xVec, f(xVec), '-b', 'Linewidth', 1);
[x0,y] = ginput(1);
plot(x0, f(x0), 'r*');

x = fzero(f, x0);
plot(x, f(x), 'g*');
hold off

disp(x);
%% Block 3 (task 3).
clc; clear; format compact; close all

%-----Parameters-------
f = @(x) x.*cos(log(x));
n = 10000;
a = 0.001;
b = 3000;
%----End of parameters-----

xVec = linspace(a, b, n); 
rootVec = zeros(1, n);
for i = 1:n
    [rootVec(i),fval, exitflag, output] = fzero(f, xVec(i));
    if(isnan(rootVec(i)))
        i
    end
end

fig_id = figure;
fig_id.Position = [300 100 1000 650];
hold on
grid on
plot(xVec, rootVec, '-r', 'Linewidth', 1);
plot(xVec, f(xVec), '-b', 'Linewidth', 1);
hold off


%% Block 4 (task 4.1)
clc; clear; format compact; close all

%-----Parameters----------
n = 1000;
mu = 5;
sigma = 2;
f = @(x) 1/(sqrt(2*pi)) .* exp(-1/2 .* x.^2);
%-----End of parameters---

x = 1:n;
y = sigma * randn(1,n) + mu;
avy = mean(y);
y1 = y * triu(ones(n))./ x;

figure;
grid on;
hold on;
plot(x, avy * ones(1,n), 'b-', x, y1, 'r-');
legend({'$\overline{X}$', '$\frac{S_n}{n}$'}, "Interpreter", "latex", 'FontSize', 15, 'Location', 'northeast');
hold off;

m = -5:0.5:5;
y2 = zeros(1,length(m));
aVec = zeros(1,n);

for i = 1:n
    y = sigma * randn(1,n) + mu;
    y1 = sum(y)/n;
    a = 1/sigma * sqrt(n) * (y1 - mu);
    aVec(i) = a;
    j = find(a < m, 1);
    y2(j) = y2(j) + 1;
    if j ~= 1
        y2(j-1) = y2(j-1) + 1;
    end
end
y2 = y2 ./ n;
x = -5:0.01:5;
figure;
grid on;
hold on;
plot(x, f(x), 'b-', m, y2, 'r-');
legend({'$~N(0,1)$', '$\frac{\sqrt{n}}{\sigma}\left( \frac{S_n}{n}-\mu \right)$'}, "Interpreter", "latex", 'FontSize', 15, 'Location', 'northeast');
hold off;
%% (task 4.2)
clc; clear; format compact; close all

%-----Parameters----------
a = 5;
b = 2;
n = 1000000;
%-----End of parameters---
x = 1:n;
y = a + b * tan(pi*(rand(1, n) - 0.5));
y1 = cumsum(y)./ x;
avy = mean(y);
figure;
grid on;
hold on;
plot(x, avy * ones(1,n), 'b-', x, y1, 'r-');
legend({'$\overline{X}$', '$\frac{S_n}{n}$'}, "Interpreter", "latex", 'FontSize', 15, 'Location', 'northeast');
hold off;

%% Block 5 (task 5)
clc; clear; format compact; close all
%-----Parameters----------
beta1 = 0.001;
beta2 = 0.0003;
alpha = 0.2;
gamma = 0.05;
delta = 0.5;
k = 0.3;
T = 1000;
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
x0Vec = [100,100,15,10];
%-----End of parameters---

sys = @(t,x) funcSIR(t, x, beta1, beta2, alpha, gamma, delta, k);

[timeVec, xMat] = ode45(sys, [0 T], x0Vec, options);
figure;
hold on;
grid on;
plot(timeVec, xMat(:,1),'b-');
plot(timeVec, xMat(:,2),'g-');
plot(timeVec, xMat(:,3),'r-');
plot(timeVec, xMat(:,4),'k-');
plot(timeVec, delta.*xMat(:,4),'m-');
legend({'Не привит', 'Привит', 'Инфицирован', 'Изолирован', 'Умер'}, 'FontSize', 15, 'Location', 'northeast');
hold off;


%% Block 6 (task 6)
clc; clear; format compact; close all

%-----Parameters----------
n = 100;
ay = -10;
by = 20;
ax = -10;
bx = 20;
alpha = 0.25;
t0 = 0;
T = 100;
h = 0.5;
x0Vec = [0, 15, 0, 12];
%-----End of parameters---

xSquare = linspace(ax,bx,n);
ySquare = linspace(ay,by,n);

fEvent = @(t,x) funcEvent6(t,x,ax,bx,ay,by);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Events', fEvent);
sys = @(t,x) func6(t, x, alpha);
t = t0;
T1 = 0;
flag = 0;
i = 0;
Video1 = struct('cdata',[],'colormap',[]);
while t <= T
    i = i + 1;
    [timeVec, xMat, te, xe, ie] = ode45(sys, [t, t+h], x0Vec, options);
    if isempty(te)
        t = t + h;
        x0Vec = xMat(end,:);
        flag = length(xMat);
    else
        t = te(1);
        x0Vec = xe;
        if ((abs(x0Vec(1)-ax)<0.001)||(abs(x0Vec(1)-bx)<0.001))&&((abs(x0Vec(3)-ay)<0.001)||(abs(x0Vec(3)-by)<0.001))
            x0Vec(2) = -x0Vec(2);
            x0Vec(4) = -x0Vec(4);
        end
        if (abs(x0Vec(1)-ax)<0.001)
            x0Vec(2) = -x0Vec(2);
            x0Vec(1) = x0Vec(1)+0.1;
        end
        if (abs(x0Vec(1)-bx)<0.001)
            x0Vec(2) = -x0Vec(2);
            x0Vec(1) = x0Vec(1)-0.1;
        end
        if (abs(x0Vec(3)-ay)<0.001)
            x0Vec(4) = -x0Vec(4);
            x0Vec(3) = x0Vec(3)+0.1;
        end
        if (abs(x0Vec(3)-by)<0.001)
            x0Vec(4) = -x0Vec(4);
            x0Vec(3) = x0Vec(3)-0.1;
        end
        flag = find(abs(timeVec-te(1))<eps, 1);
    end
    
    p_id = plot(xMat(1:flag,1), xMat(1:flag,3), 'b-', xMat(flag,1), xMat(flag,3), 'ro', xSquare, ones(1,n)*ay, 'k-', xSquare, ones(1,n)*by, 'k-',...
    ones(1,n)*ax, ySquare, 'k-', ones(1,n)*bx, ySquare, 'k-');
    axis equal;
    set(gca,'FontSize',14);
    hold off;
    Video1(i) = getframe(gcf);
    delete(p_id);
end


vidObj = VideoWriter('Video36.avi');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);
writeVideo(vidObj,Video1);
close(vidObj);


%% Block 7 (task 7)
clc; clear; format compact; close all

%-----Parameters----------

%спираль
% m1 = 10;
%красный
% m2 = 1000;
% G = 1;
% T = 100;
% t0 = 0;
% h = 0.5;
% x0Vec = [0, 0, 0, 10, 10, -1, 0, -2];

%вокруг общего центра
% m1 = 10;
% m2 = 10;
% G = 1;
% T = 100;
% t0 = 0;
% h = 0.5;
% x0Vec = [-3.5, 0, 0, 1, 3.5, 0, 0, -1];

% m1 = 120;
% m2 = 100;
% G = 0.367;
% T = 100;
% t0 = 0;
% h = 0.1;
% x0Vec = [0, 5, 0, 5, 0, -5, 0.4, -5];

m1 = 120;
m2 = 100;
G = 0.367;
T = 100;
t0 = 0;
h = 0.1;
x0Vec = [0, 5, 0, 5, 0, -5, 0.4, -5];

%-----End of parameters---

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sys = @(t,x) func7(t, x, m1, m2, G);
t = t0;
T1 = 0;
i = 0;
Video1 = struct('cdata',[],'colormap',[]);
while t <= T
    i = i + 1;
    [timeVec, xMat] = ode45(sys, [t, t+h], x0Vec, options);
    t = t + h;
    x0Vec = xMat(end,:);
    hold on;
    p_id = plot(xMat(:,1), xMat(:,3), 'b-', xMat(:,5), xMat(:,7), 'r-');
    axis equal;
    set(gca,'FontSize',14);
    hold off;
    Video1(i) = getframe(gcf);
    pause(0.1);
%     delete(p_id);
end


vidObj = VideoWriter('Video37.avi');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);
writeVideo(vidObj,Video1);
close(vidObj);

%% Block 8 (task 8)
clc; clear; format compact; close all

%-----Parameters----------
a = -10;
b = 10;
%седло
%aVec = [1,1,1,-1];

%фокус
%aVec = [-1,1,-1,-1];

%центр
%aVec = [-1,1,-2,1];

%узел
%aVec = [2,-1,-1,2];

%-----End of parameters---

[X,Y] = meshgrid(a:b,a:b);
U = aVec(1)*X + aVec(2)*Y;
V = aVec(3)*X + aVec(4)*Y;
hold on;
grid on;
quiver(X,Y,U./sqrt(U.^2+V.^2),V./sqrt(U.^2+V.^2));
hold off;

%% Block 9 (task 9)
clc; clear; format compact; close all

%-----Parameters----------
n = 5;
m = 5;
%-----End of parameters---
P = GenerateTable(n,m);
%%
DesperateSapper(P);

%% Block 10 (task 10)
clc; clear; format compact; close all

%-----Parameters----------
t0 = 0;
T = pi;
options = bvpset('RelTol',1e-6,'AbsTol',1e-6);
%-----End of parameters---

f = @(t,x) func10(t,x);
bc = @(xl,xr) func10V(xl,xr);
t_meshVec = linspace(t0,T,100);
solinit = bvpinit(t_meshVec,[1,0]);

sol = bvp5c(f,bc,solinit,options);
timeVec = sol.x;
xMat = sol.y;

figure;
hold on;
grid on;
plot(timeVec, xMat(1,:),'b-');
plot(timeVec, xMat(2,:),'r-');
xlabel('time');
ylabel('x(t)');
legend({'$x_1$','$x_2$'},'Interpreter','latex');
set(gca,'Fontsize',16);
hold off;


%% Block 11 (task 11)
clc; clear; format compact; close all
%-----Parameters----------
a = 0;
b = 1;
%-----End of parameters---

