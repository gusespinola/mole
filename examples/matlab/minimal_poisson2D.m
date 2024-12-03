% 2D Staggering example using a 2D Mimetic laplacian

clc
close all

addpath('../../src/matlab')

k = 2; % Order of accuracy
%h = 0.02; % Discretization size
%m = 1/h; % Vertical resolution
%n = 1/h; % Horizontal resolution
m = 32;
n = 32;
h = 1/m;

tic
L = lap2D(k, m, 1, n, 1); % 2D Mimetic laplacian operator
L2 = L + robinBC2D(k, m, 1, n, 1, 1, 0); % Dirichlet BC

RHS = zeros(m+2, n+2);

RHS(1, :) = 100; % Known value at the bottom boundary

RHS = reshape(RHS, [], 1);

SOL = L2\RHS;

SOL = reshape(SOL, m+2, n+2);
toc
figure(1)
imagesc(SOL)
title('2D Poisson''s equation - numerical solution')
xlabel('M_x')
ylabel('M_y')
set(gca, 'YDir', 'Normal')
colorbar
hold on

% Analytical solution - Fourier series
x=0:h:1;
y=0:h:1;
[X, Y]=meshgrid(x, y);
uxy=zeros(size(X));
tic
for n=1:100
    cn=(1-(-1)^n)/n;
    uxy=uxy+cn*sin(n*pi*X)...
        .*((exp(n*pi*Y)-exp(n*pi*(2.-Y))))...
        /(1-exp(2*n*pi));
end
uxy=(200/pi).*uxy;
toc
%surf(X, Y, uxy)
figure(2)
imagesc(flipud(uxy))
title('2D Poisson''s equation - analytical solution')
xlabel('M_x')
ylabel('M_y')
colorbar
set(gca, 'CLim', [0 100])
hold on