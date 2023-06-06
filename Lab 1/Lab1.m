%% Block 1 (task 1).
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
plot(xVec, fVec, 'b-', 'Linewidth', 1);
xlabel('$x$', 'Interpreter','latex');
ylabel('$f(x)$', 'Interpreter','latex');
xlim([a - 0.5, b + 0.5]);
ylim([min(fVec) - 1, max(fVec) + 1]);

[~, I_max] = find(abs(fVec - max(fVec)) < 0.0001);
[~, I_min] = find(abs(fVec - min(fVec)) < 0.0001);
plot(xVec(I_max), fVec(I_max), '.', 'Markersize', 15, 'Color', [0, 0, 1], 'Linewidth', 1);
plot(xVec(I_min), fVec(I_min), '.', 'Markersize', 15, 'Color', [1, 0, 0], 'Linewidth', 1);

hold off

legend({'$f(x)$', '$max$', '$min$'}, "Interpreter", "latex", 'FontSize', 15, 'Location', 'northeast');

%% Block 2 (task 2).
clc; clear; format compact; close all

%-----Parameters-------
n = 4;
%----End of parameters-----
if ((n == round(n)) && (n > 0))

    Vec = 9:18:n;
    disp(Vec);
    
    
    a = 1:n;
    A_Mat = repmat(a', 1, n);
    disp(A_Mat);
    
    b = 1:(n*(n+1));
    B_Mat = reshape(b, n+1, n);
    disp(B_Mat');

   % c = reshape(B_Mat', n*(n+1), 1);
    c = B_Mat(1:2:end);
    disp(c);
    
    
end

%% Block 3 (task 3).
clc; clear; format compact; close all

%-----Parameters-------
n = 4;
%----End of parameters-----

X_Mat = round(ones(n)+99*rand(n));

disp(X_Mat);
disp(max(diag(X_Mat)));
Y_Mat = sortrows(X_Mat, 3);
disp(Y_Mat);

[~, order] = sort(X_Mat(:,3));
Sort_Mat = X_Mat(order,:);
disp(Sort_Mat);

%% Block 4 (task 4).
clc; clear; format compact; close all

%-----Parameters-------
n = 4;
m = 5;
%----End of parameters-----

x = ceil(10*rand(n, 1));
y = ceil(10*rand(1, m));
Mat = x*y;
disp(Mat);


%% Block 5 (task 5).
clc; clear; format compact; close all

%-----Parameters-------
n = 19;
%----End of parameters-----

if isprime(n) == 1
    disp('prime');
else
    disp('not prime');
end


A = rand(n);
b = rand(n,1);
if abs(det(A)) > 0.0001
    x1 = A^(-1)*b;
    x2 = A^(-1)*b;
    x3 =  A\b;
    [L, U] = lu(A);
    y = L\b;
    x4 = U\y;
    ans1 = norm(x1 - x2);
    ans2 = norm(x2 - x3);
    ans3 = norm(x1 - x3);
    ans4 = norm(x1 - x4);
    ans5 = norm(x2 - x4);
    ans6 = norm(x3 - x4);
else
    disp('matrix is degenerate!')
end



%% Block 6 (task 6).
clc; clear; format compact; close all

%-----Parameters-------
n = 4;
m = 5;
%----End of parameters-----

a = ceil(10*rand(n, 1));
b = ceil(10*rand(m, 1));
res = max(max(a) - min(b), max(b) - min(a));
disp(res);



%% Block 7 (task 7).
clc; clear; format compact; close all

%-----Parameters-------
n = 4;
k = 3;
%----End of parameters-----

A = ceil(5*rand(n, k));
temp = A*A';
dtemp = diag(temp);
rMat = repmat(dtemp, 1, n);
distant = sqrt(rMat - 2*temp + rMat');
fprintf('res = \n')
disp(distant);



%% Block 8 (task 8).
clc; clear; format compact; close all

%-----Parameters-------
n = 4;
m = 5;
%----End of parameters-----

X_Mat = ceil(100*rand(m, n));
Y_Mat = diag(diag(X_Mat));
y = reshape(Y_Mat, [], 1);
id = find(y == 0);
y(id) = 5;
Y_Mat = reshape(y, min(n, m), min(n, m));
if n == min(n,m)
    Y_Mat = [Y_Mat;  5*ones(m-n, n)];
else
    Y_Mat = [Y_Mat 5*ones(m, n-m)];
end

disp(Y_Mat);



%% Block 9 (task 9).
clc; clear; format compact; close all

%-----Parameters-------
n = 4;
%----End of parameters-----

A = ceil(10*rand(n));
B = ceil(10*rand(n));

nVec = 1:100:10000;
TVec1 = time1(nVec);
TVec2 = time2(nVec);

makeaplot(nVec, TVec1);
makeaplot(nVec , TVec2);


%% Block 10 (task 10).
clc; clear; format compact; close all
  
X = [NaN 1 2; NaN 0 6; 1 5 NaN];

Y = nanMean(X);
disp(Y);


%% Block 11 (task 11).
clc; clear; format compact; close all

%-----Parameters-------
n = 10000000;
a = 5;
sigma = 2;

%----End of parameters-----

vec = sigma*randn(n, 1) + a;

prop = length(find((vec <= (a + 3*sigma))&(vec >= a - 3*sigma)))./n;


%% Block 12 (task 12).
clc; clear; format compact; close all

%-----Parameters-------
a = -100;
b = 100;
n = 1000;
h = 0.001;
f = @(x) sin(x)./x;
%----End of parameters-----

X = (a:h:b)';
Y = (f(X))';
Y(isnan(Y)) = 1;

tic
int1 = trapz(X, Y)
t = toc
tic
int2 = simpson(X, f)
t = toc
tic
int3 = rectangles(X, f)
t = toc


X1 = linspace(a, b, n);

antider1 = cumtrapz(X1, f(X1));
antider2 = simpsonn(X1, f);
antider3 = rectangless(X1, f);


fig_id = figure;
fig_id.Position = [300 100 900 550];
hold on
grid on
title('Integral')
plot(X1,antider1,'r-')
plot(X1,antider2,'b-')
plot(X1,antider3,'go')
legend('trapz', 'simpson', 'rectangles')
xlabel('$x$','Interpreter','Latex')
ylabel('$integral$','Interpreter','Latex')
hold off




h = 0.0001:0.0001:0.1;
err1 = zeros(length(h));
err2 = zeros(length(h));
err3 = zeros(length(h));

for i = 1:length(h)
    X = a:h(i):b;
    X1 = a:(h(i)/2):b;
    Y = f(X);
    Y1 = f(X1);
    Y(isnan(Y)) = 1;
    err1(i) = abs(simpson(X, f) - simpson(X1, f));
    err2(i) = abs(rectangles(X, f) - rectangles(X1, f));
    err3(i) = abs(trapz(X, Y) - trapz(X1, Y1));
end

fig_id = figure;
fig_id.Position = [300 100 1000 650];

subplot(2, 2, 1)
hold on
grid on
plot(h, err1, 'b-')
hold off

subplot(2, 2, 2)
hold on
grid on
plot(h, err2, 'g-')
hold off

subplot(2, 2, 3)
hold on
grid on
plot(h, err3, 'r-')
hold off

%% Block 13 (task 13).
clc; clear; format compact; close all

%-----Parameters-------
a = -3*pi;
b = 3*pi;
n = 12000;

%----End of parameters-----

f = @(x) sin(x);
fder = @(x) cos(x);
h = logspace(-10,-2);
x0 = 14;
fder_c = (f(x0+h)-f(x0-h))./(2.*h);
fder_r = (f(x0+h)-f(x0))./(h);


Y1Vec = abs(fder(x0)-fder_c);
Y2Vec = abs(fder(x0)-fder_r);


fig_id = figure;
fig_id.Position = [300 100 1000 650];

hold on

grid on
loglog(h, Y1Vec, 'b-')
loglog(h, Y2Vec, 'r-')
set(gca, 'XScale', 'log', 'YScale', 'log')
hold off





