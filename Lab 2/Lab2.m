  %% Block 1 (task 1).
clc; clear; format compact; close all

%-----Parameters-------
a = -1;
b = 1;
n = 1000;
m = 5000;
%f = @(x) 1./(1+25*x.^2);
%f = @(x) sign(x);
f = @(x) sin(x)./x;
%----End of parameters-----


x = linspace(a, b, n);
xx = linspace(a, b, m);
compareInterp(x,xx,f);



%% Block 2 (task 2).
clc; clear; format compact; close all

x = linspace(0, 1, 1000);

%M
xMVec = [x zeros(1, 1000) x];
yMVec = [(1-x) x zeros(1, 1000)];

%N
xNVec = [x+2 3*ones(1, 1000) (x+2)];
yNVec = [x x zeros(1, 1000)];

Video1 = struct('cdata', [], 'colormap', []);
%1
flag = 1;
rVec = 0:0.05:10;
i = 0;
fig_id = figure;
fig_id.Position = [300 100 1000 650];
while flag
    t = 100;
    xrMat = repmat(xMVec, t, 1);
    yrMat = repmat(yMVec, t, 1);
    i = i+1;
    r = rVec(i);
    theta = linspace(0, 2*pi, t);
    xr = r*cos(theta);
    yr = r*sin(theta);
    xrMat = xrMat + repmat(xr', 1, length(xMVec));
    yrMat = yrMat + repmat(yr', 1, length(yMVec));

    posVec = convhull(xrMat, yrMat);
    xVec = xrMat(posVec);
    yVec = yrMat(posVec);

    hold on
    grid on   

    p_id = plot(xVec, yVec,  xNVec, yNVec);
    if (inpolygon(xNVec, yNVec, xVec, yVec))
        flag = 0;
    end
    str = ['r = ', num2str(r)];
    title(str)
    set(gca, 'Fontsize', 14)
    hold off

    Video1(i) = getframe(gcf);
    delete(p_id)
end
rmn = r;

%2
flag = 1;
rVec = 0:0.05:10;
j = 0;
while flag
    t = 100;
    xrMat = repmat(xNVec, t, 1);
    yrMat = repmat(yNVec, t, 1);
    j = j+1;
    r = rVec(j);
    theta = linspace(0, 2*pi, t);
    xr = r*cos(theta);
    yr = r*sin(theta);
    xrMat = xrMat + repmat(xr', 1, length(xNVec));
    yrMat = yrMat + repmat(yr', 1, length(yNVec));

    posVec = convhull(xrMat, yrMat);
    xVec = xrMat(posVec);
    yVec = yrMat(posVec);
    
   
    hold on
    grid on   

    p_id = plot(xVec, yVec,  xMVec, yMVec);
    if (inpolygon(xMVec, yMVec, xVec, yVec))
        flag = 0;
    end
    str = ['r = ', num2str(r)];
    title(str)
    set(gca, 'Fontsize', 14)
    hold off
 

    Video1(j+i) = getframe(gcf);
    delete(p_id)
end

rnm = r;


vid0bj = VideoWriter('myVideo.avi');
vid0bj.FrameRate = 10;
vid0bj.Quality = 100;
open(vid0bj);
writeVideo(vid0bj, Video1);
close(vid0bj);






%% Block 3 (task 3).
clc; clear; format compact; close all

%-----Parameters-------
a = 0;
b = 5;
n = 1000;
%----End of parameters-----

%fn = @(n, x) n.*x./(n+x.^n);
%f = @(x) x;

fn = @(n, x) n.*x.^3./(1+n.^2.*x.^6);
f = @(x) 0.*x;


pl = convergenceFunc(fn,f,a,b,n);

%% Block 4 (task 4).
clc; clear; format compact; close all

%-----Parameters-------
f = @(x) round(5*x);
n = 50;
%----End of parameters-----

chebApprox(f, n);

%% Block 5 (task 5).
clc; clear; format compact; close all

%-----Parameters-------
a = -pi;
b = pi;
n = 1000;
%----End of parameters-----

f = @(x) sin(log(1 + abs(x))- x.^2);

xVec = linspace(a, b, n);
fVec = f(xVec);


fig_id = figure;
fig_id.Position = [300 100 1000 650];

hold on

grid on
title('$f(x) = \sin(\ln (1+|x|)-x^2)$', 'Interpreter', 'latex');
plot(xVec, fVec, 'g-', 'Linewidth', 1);
xlabel('$x$', 'Interpreter','latex');
ylabel('$f(x)$', 'Interpreter','latex');
xlim([a - 0.5, b + 0.5]);
ylim([min(fVec) - 1, max(fVec) + 1]);

[f_max, I_max] = max(fVec);
d = diff(fVec);
d(d<0) = -1;
d(d>0) = 1;
d = diff(d);
I_min_0 = find(d > 0);
[~, temp] =  min((xVec(I_max) - xVec(I_min_0)).^2 + (fVec(I_max)-fVec(I_min_0)).^2);
I_min = I_min_0(temp);

plot(xVec(I_max), fVec(I_max), '.', 'Markersize', 15, 'Color', [0, 0, 1], 'Linewidth', 1);
plot(xVec(I_min_0), fVec(I_min_0), '.', 'Markersize', 15, 'Color', [1, 0, 0], 'Linewidth', 1);
if (I_max>I_min)
    h = -1;
else
    h = 1;
end
comet(xVec(I_max:h:I_min), fVec(I_max:h:I_min), 0.9);
hold off

legend({'$f(x)$', '$max$', '$min$'}, "Interpreter", "latex", 'FontSize', 15, 'Location', 'northeast');
%% Block 6 (task 6).
clc; clear; format compact; close all

%-----Parameters-------
h = 1;
r = 1;
x = @(t) r*t - h*sin(t);
y = @(t) r - h*cos(t);
t0 = 0;
t1 = 2*pi;
N = 6;
%----End of parameters-----
getEqual(x,y,t0,t1,N)
%% Block 7 (task 7).
clc; clear; format compact; close all

%-----Parameters-------
x0 = 0;
y0 = 1;
h = 1.5;
r = 1;
T = 30;
n = 100;
%----End of parameters-----

Trochoid(x0,r,r,h,T,n);

%% Block 8 (task 8).
clc; clear; format compact; close all


data = readtable('2020_2.csv');
IDVec = data.ID;
LVec = data.CREDPAY;
terVec = data.MEST;
hhamtVec = data.CHLICO;
SCell = data.DOXODSN;
CCell = data.DENRAS;
size = length(IDVec);
CVec = zeros(size, 1);
SVec = zeros(size, 1);


for i = 1:size
    c = CCell{i};
    s = SCell{i};
    c(c == ',') = '.';
    s(s == ',') = '.';
    CVec(i) = str2double(c);
    SVec(i) = str2double(s);
end

AMat = [IDVec,  LVec,  terVec, hhamtVec,  SVec, CVec];
[str, ~] = find(isnan(AMat));
AMat = AMat(~ismember(1:size, str), :);

%% Block 8 (task 8).
clc; clear; format compact; close all
%-----Parameters----------
a = -1;
b = 3;
n = 100000;
%-----End of parameters---

fid = fopen('2020_2.csv');
tline = fgetl(fid)
q = 0;
aMat = zeros(n,6);
while ~feof(fid)
    tline = fgetl(fid);
    k = find(tline==';');
    m = find(tline==',');
    id = str2num(tline(1:k(1)));
    st = tline(k(1)+1:k(2)-1);
    k1 = find(st==',');
    if k1>0
        st(k1)='.';
    end
    credpay = str2double(st);
    mest = str2num(tline(k(2):k(3)));
    chlico = str2num(tline(k(3):k(4)));
    st = tline(k(4)+1:k(5)-1);
    k1 = find(st==',');
    if k1>0
        st(k1)='.';
    end
    doxod = str2double(st);
    st = tline(k(5)+1:length(tline));
    k1 = find(st==',');
    if k1>0
        st(k1) = '.';
    end
    denras = str2double(st);
    
    if (doxod == 0)||(denras == 0)||(chlico>15)||(doxod/(3.*chlico) > 50000)||(denras/(3*chlico) > 50000)
        continue;
    elseif (isnan(id))||(isnan(credpay))||(isnan(mest))||(isnan(chlico))||(isnan(doxod))||(isnan(denras))
        continue;
    end
    q = q+1;
    aMat(q,1) = id;
    aMat(q,2) = credpay;
    aMat(q,3) = mest;
    aMat(q,4) = chlico;
    aMat(q,5) = doxod;
    aMat(q,6) = denras;
  
end
fclose(fid);
av1 = find(aMat(1:q,3)==1);
av2 = find(aMat(1:q,3)==2);

bMat = zeros(4,2);
bMat(1,1) = mean(aMat(av1,5)./aMat(av1,4));
bMat(2,1) = median(aMat(av1,5)./aMat(av1,4));
bMat(3,1) = mean(aMat(av1,6)./aMat(av1,4));
bMat(4,1) = median(aMat(av1,6)./aMat(av1,4));

bMat(1,2) = mean(aMat(av2,5)./aMat(av2,4));
bMat(2,2) = median(aMat(av2,5)./aMat(av2,4));
bMat(3,2) = mean(aMat(av2,6)./aMat(av2,4));
bMat(4,2) = median(aMat(av2,6)./aMat(av2,4));


fid = fopen('inc_cons.txt','w');
fprintf(fid,'Средний доход в городе: %.2f\n', bMat(1,1));
fprintf(fid,'Медианный доход в городе: %.2f\n', bMat(2,1));
fprintf(fid,'Средний расход в городе: %.2f\n', bMat(3,1));
fprintf(fid,'Медианный расход в городе: %.2f\n\n', bMat(4,1));

fprintf(fid,'Средний доход в селе: %.2f\n', bMat(1,2));
fprintf(fid,'Медианный доход в селе: %.2f\n', bMat(2,2));
fprintf(fid,'Средний расход в селе: %.2f\n', bMat(3,2));
fprintf(fid,'Медианный расход в селе: %.2f\n', bMat(4,2));

fclose(fid);
%%
xVec = categorical({'av\_d','med\_d','av\_r','med\_r'});
figure;
hold on;
grid on;
b = bar(xVec,bMat./3);
b(1).FaceColor = 'Flat';
b(2).FaceColor = 'Flat';
min_pos = find(abs(min(min(bMat))-bMat)<eps);
max_pos = find(abs(max(max(bMat))-bMat)<eps);

ch1 = floor(min_pos/4);
ch2 = mod(min_pos,4);
if (ch1==0)
    ch1 = 1;
end
if (ch2==0)
    ch2 = 4;
end
b(ch1).CData(ch2,:) = [1,1,0];
ch1 = floor(max_pos/4);
ch2 = mod(max_pos,4);
if (ch1==0)
    ch1 = 1;
end
if (ch2==0)
    ch2 = 4;
end
b(ch1).CData(ch2,:) = [0,1,0];

set(gca,'FontSize',16);
hold off;

zg = length(find((aMat(1:q,2)>0)&(aMat(1:q,3)==1)));
zs = length(find((aMat(1:q,2)>0)&(aMat(1:q,3)==2)));

figure;
subplot(1,2,1);
pie([zg, length(av1)-zg]);
title('Город');

subplot(1,2,2);
pie([zs, length(av2)-zs]);
title('Село');

zaem1 = find(aMat(av1,2)>0);
nezaem1 = find(aMat(av1,2)<eps);
zaem2 = find(aMat(av2,2)>0);
nezaem2 = find(aMat(av2,2)<eps);

fid = fopen('table.txt','w');
fprintf(fid,'Средний доход заемщиков в городе: %.2f\n', mean(aMat(zaem1,5)./aMat(zaem1,4))/3);
fprintf(fid,'Средний расход заемщиков в городе: %.2f\n', mean(aMat(zaem1,6)./aMat(zaem1,4))/3);
fprintf(fid,'Средний состав семьи заемщиков в городе: %.2f\n\n', mean(aMat(zaem1,4)));

fprintf(fid,'Средний доход остальных в городе: %.2f\n', mean(aMat(nezaem1,5)./aMat(nezaem1,4))/3);
fprintf(fid,'Средний расход остальных в городе: %.2f\n', mean(aMat(nezaem1,6)./aMat(nezaem1,4))/3);
fprintf(fid,'Средний состав семьи остальных в городе: %.2f\n\n', mean(aMat(nezaem1,4)));

fprintf(fid,'Средний доход заемщиков в селе: %.2f\n', mean(aMat(zaem2,5)./aMat(zaem2,4))/3);
fprintf(fid,'Средний расход заемщиков в селе: %.2f\n', mean(aMat(zaem2,6)./aMat(zaem2,4))/3);
fprintf(fid,'Средний состав семьи заемщиков в селе: %.2f\n\n', mean(aMat(zaem2,4)));

fprintf(fid,'Средний доход остальных в селе: %.2f\n', mean(aMat(nezaem2,5)./aMat(nezaem2,4))/3);
fprintf(fid,'Средний расход остальных в селе: %.2f\n', mean(aMat(nezaem2,6)./aMat(nezaem2,4))/3);
fprintf(fid,'Средний состав семьи остальных в селе: %.2f\n\n', mean(aMat(nezaem2,4)));

fclose(fid);

%% Block 9 (task 9)
clc
format compact
aMat = aMat(1:q,1:6);
min(aMat(:,5)./aMat(:,4))/10000
max(aMat(:,5)./aMat(:,4))/10000
min(aMat(:,6)./aMat(:,4))/10000
max(aMat(:,6)./aMat(:,4))/10000


[xMat,yMat] = meshgrid(0:1000:0000,0:1000:60000);
fMat = zeros(size(xMat));

% for m = 0:1000:1000000
%     i = i+1;
%     a = length(find((aMat(:,5)./aMat(:,4)>=m)&(aMat(:,5)./aMat(:,4)<=m+1000)&(aMat(:,3)==1)));
%     b = length(find((aMat(:,6)./aMat(:,4)>=m)&(aMat(:,6)./aMat(:,4)<=m+1000)&(aMat(:,3)==1)));
%     xVec(i) = a;
%     yVec(i) = b;
% end

for i=1:q
    if aMat(i,3)==1
        x1 = floor(aMat(i,5)./aMat(i,4)/1000)+1
        x2 = ceil(aMat(i,5)./aMat(i,4)/1000)+1;
        y1 = floor(aMat(i,6)./aMat(i,4)/1000)+1;
        y2 = ceil(aMat(i,6)./aMat(i,4)/1000)+1;
        fMat(x1,y1) = fMat(x1,y1)+1;
        fMat(x1,y2) = fMat(x1,y2)+1;
        fMat(x2,y1) = fMat(x2,y1)+1;
        fMat(x2,y2) = fMat(x2,y2)+1;
    end
end

zMat = smoothn(fMat);
xVec = 0:1000:50000;
yVec=xVec;
I = trapz(yVec,trapz(xVec,zMat,2));


figure;
hold on;
grid on;
surf(xMat,yMat,zMat/I, 'EdgeColor', 'none');
xlabel('$DOXOD$','Interpreter','latex','FontSize', 14);
ylabel('$RASXOD$','Interpreter','latex','FontSize', 14);
zlabel('$NUM$','Interpreter','latex','FontSize', 14);
title('GOROD');
view(3);
hold off;

for i=1:q
    if aMat(i,3)==2
        x1 = floor(aMat(i,5)./aMat(i,4)/1000)+1;
        x2 = ceil(aMat(i,5)./aMat(i,4)/1000)+1;
        y1 = floor(aMat(i,6)./aMat(i,4)/1000)+1;
        y2 = ceil(aMat(i,6)./aMat(i,4)/1000)+1;
        fMat(x1,y1) = fMat(x1,y1)+1;
        fMat(x1,y2) = fMat(x1,y2)+1;
        fMat(x2,y1) = fMat(x2,y1)+1;
        fMat(x2,y2) = fMat(x2,y2)+1;
    end
end
figure;
hold on;
grid on;
surf(xMat,yMat,smoothn(fMat)/I, 'EdgeColor', 'none');
xlabel('$DOXOD$','Interpreter','latex','FontSize', 14);
ylabel('$RASXOD$','Interpreter','latex','FontSize', 14);
zlabel('$NUM$','Interpreter','latex','FontSize', 14);
title('SELO');
view(3);
hold off;


%% Block 10 (task 10)
clc; clear; format compact; close all

%-----Parameters----------
a = 0;
b = 0;
q = 100;
n = 100;
%-----End of parameters---

xVec = linspace(-2,5,q);
yVec = linspace(-2,5,q);
[xMat, yMat] = meshgrid(xVec,yVec);
f = @(x,y,a,b) sin(x+a)+cos(y+b);

Video1(1:n) = struct('cdata',[],'colormap',[]);

for i = 1:n
    a = a + pi/4;
    b = b + pi/6;
    fMat = f(xMat,yMat,a,b);
    posmin = pov_min(fMat);
    posmax = pov_max(fMat);
    hold on;
    grid on;
    s_id = surf(xMat,yMat, fMat,'EdgeColor','none');
    sc1_id = scatter3(xMat(posmin),yMat(posmin),fMat(posmin),'MarkerFaceColor',[1 0 0]);
    sc2_id = scatter3(xMat(posmax),yMat(posmax),fMat(posmax),'MarkerFaceColor',[0 0 1]);
    view(3);
    set(gca,'FontSize',14);
    hold off;
            
    Video1(i) = getframe(gcf);
    if i~=n
        delete(s_id);
        delete(sc1_id);
        delete(sc2_id);
    end
end
save('Test.mat','Video1');
%%
clear
load('Test.mat','Video1');

movie(Video1);
%%
hold on;
grid on;
contour(xMat,yMat,fMat);
colorbar;
set(gca, 'FontSize', 16);
hold off;
%%
vidObj = VideoWriter('Video10.mat');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);
writeVideo(vidObj,Video1);
close(vidObj);
%%
vidObj = VideoWriter('Video10.avi');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);
writeVideo(vidObj,Video1);
close(vidObj);

%% Block 12 (task 12)
clc; clear; format compact; close all
%-----Parameters----------
R = 400;
N = 5;
%-----End of parameters---

%shaperead, geoshow, distance, km2deg, deg2km, plotm

ax = worldmap('madagascar');
land = shaperead('landareas.shp', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5]);

%%
colormap gray

[lat, lon] = inputm(N);

t = 0:0.01:2*pi;
rxMat = repmat(lat, 1, length(t));
ryMat = repmat(lon, 1, length(t));
rxMat = rxMat + repmat(km2deg(L*cos(t), 'madagascar'), N, 1);
ryMat = ryMat + repmat(km2deg(L*sin(t), 'madagascar'), N, 1);

plotm(rxMat',ryMat','r-');

arclen = distance([lat(1:end-1),lon(1:end-1)],[lat(2:end),lon(2:end)]);
rlen = deg2km(arclen, 'madagascar')
if ~isempty(find(rlen>2*L, 1))
    disp('NO');
else
    disp('YES');
end

%% Block 13 (task 13)
clc; clear; format compact; close all

%-----Parameters----------
alpha = inf;
params.num = 100;
params.a = -10;
params.b = 10;
params.lev = 5;
params.color = 'blue';
%-----End of parametrs---

drawSphere(alpha, params);

%% Block 14 (task 14)
clc; clear; format compact; close all

%-----Parameters----------
alphas = [0.5,1.5,2,inf];
colors = ["black","red","blue","red"];
edges = ["none","none","none","none"];
trans = [1,0.6,0.4,0.2];
%-----End of parameters---

drawManySpheres(alphas,colors,edges,trans);