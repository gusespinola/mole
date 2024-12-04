% 2D Staggering example using a 2D Mimetic laplacian

clc
close all

addpath('../../src/matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 4; % Order of accuracy
%h = 0.02; % Discretization size
%m = 1/h; % Vertical resolution
%n = 1/h; % Horizontal resolution
mx = 256;
ny = 256;
h = 1/mx;
hy = 1/ny;

L2 = lap2D(k, mx, 1, ny, 1); % 2D Mimetic laplacian operator
L = L2 + robinBC2D(k, mx, 1, ny, 1, 1, 0); % Dirichlet BC

% Condition number (estimated) for LHS, no preconditioning
CondLNoPrec = condest(L);

RHS = zeros(mx+2, ny+2);

% Source term
RHS(1, :) = 100; % Known value at the bottom boundary

RHS = reshape(RHS, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL SOLUTION - BACKSLASH - NO PRECONDITIONING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usual built-in MATLAB 'backslash' (\) method, no preconditioning
disp('BACKSLASH');
tic;
SOLBNP = L\RHS;
timeBNP = toc();
SOLBNPRS = reshape(SOLBNP, mx+2, ny+2);
disp(timeBNP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=0:h:1;
y=0:h:1;
[X, Y]=meshgrid(x, y);
uxy=zeros(size(X));
%tic
for n=1:100
    cn=(1-(-1)^n)/n;
    uxy=uxy+cn*sin(n*pi*X)...
        .*((exp(n*pi*Y)-exp(n*pi*(2.-Y))))...
        /(1-exp(2*n*pi));
end
uxy=(200/pi).*uxy;
%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEFT-PRECONDITIONING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('LEFT-PRECONDITIONING');
tic;
% Jacobi
M_jacobi = (diag(diag(L)))^-1;
L = M_jacobi * L;
RHS = M_jacobi * RHS;
CondLpostJacobi = condest(L);

% Gauss-Seidel
% M_GS = (tril(L1))^-1;
% L= M_GS*L1;
% RHS = M_GS*RHS1;
% CondLpostGS = condest(L);

% SOR
% n=size(L);
% tStart=tic;
% omega=1;
% D=diag(diag(L));
% L=-(tril(L)-D);
% M_sor=(1/omega)*D-L;
% L=M_sor\L;
% RHS=M_sor\RHS;

% % Preconditioner ILU for GMRES & BICG
% setup.type = 'ilutp';
% setup.droptol = 1e-6;
% setup.udiag=1;
% disp('ILU')
% tic
% [LT,UT] = ilu(sparse(L), setup);
% timeILU = toc();
% disp(timeILU)
% M=LT * UT;
% %         clear L U;
%
% disp('INVERSAS')
% tic
% % Is there an alternative?
% L=M\L;
% RHS=M\RHS;
% timeMatMul = toc();
% disp(timeMatMul)
% CondLpostILU = condest(L);

timePREC = toc();
disp(timePREC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL SOLUTION - BACKSLASH - PRECONDITIONINED SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usual built-in MATLAB 'backslash' (\) method, preconditioned system
disp('BACKSLASH');
tic;
SOLBP = L\RHS;
timeBP = toc();
SOLBPRS = reshape(SOLBP, mx+2, []);
disp(timeBP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KRYLOV SUBSPACE METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL PARAMETERS
tol = 1e-9;
maxit = 100;

% REEMPLAZAR CON cd (path_to_solvers) O EQUIVALENTE
cd('C:\Users\HP\OneDrive\MSc Thesis\krysbas-dev\test-prec-gmres\src\solvers\')

% BUILT-IN GMRES
% m = 30;
% [SOL0, flag0, relres0, iter0, resvec0] = ...
%     gmres(L, RHS, m, tol, maxit);

% BUILT-IN GMRES + JACOBI
% [SOL00, flag00, relres00, iter00, resvec00] = ...
%     gmres(L, RHS, m, tol, maxit, diag(diag(L)));

% PREC_GMRES
% disp('MI-PREC-GMRES');
% m = 30;
% print = 1;
% color = 'c-';
% flagprec = 1; 
% itermax = 1000; 
% %itermax = maxit;
% Name_Matrix = 'poisson';
% tic
% [vec11, vny1, xx1, relresPG]= mi_prec_GMRES_m(L,RHS,m,itermax,tol, color,print, Name_Matrix, flagprec);
% timePG = toc();
% disp(timePG);

% PD_GMRES
mInitial = 30;
disp('PD-GMRES');
[SOLPD, flagPD, relresvecPD, mvec, timePD] = ...
    pd_gmres(L, RHS, mInitial, [], [], tol, maxit, [], [1, 1]);
disp(timePD);
SOLPDRS = reshape(SOLPD, mx + 2, ny + 2);

% LGMRES
m = 27;
k = 3;
disp('LGMRES');
[SOLL, flagL, relresvecL, kdvecL, timeL] = ...
    lgmres(L, RHS, m, k, tol, maxit);
disp(timeL);
SOLLRS = reshape(SOLL, mx + 2, ny + 2);

% GMRES-E
m = 27;
k = 3;
disp('GMRES-E');
[SOLE, flagE, relresvecE, kdvecE, timeE] = ...
    gmres_e(L, RHS, m, k, tol, maxit);
disp(timeE);
SOLERS = reshape(SOLE, mx + 2, ny + 2);

% Parameters for built-in MATLAB iterative methods
m = 30;
maxcycles = 100; % maxcycles: max number of outer iterations
maxit = m * maxcycles; % maxit: max number of inner iterations
                        % necessary as built-in gmres() input parameter 

% GMRES(m)
disp('GMRES(m)');
tic;
[SOLG, flagG, resG, iterG, resvecG] = gmres(L, RHS, m, tol, maxcycles);
relresvecG = resvecG ./ norm(RHS);
timeG = toc();
disp(timeG);
SOLGRS = reshape(SOLG, mx + 2, ny + 2);

% GMRES(m) + JACOBI
disp('PREC-GMRES(m), JACOBI');
tic;
[SOLGP, flagGP, resGP, iterGP, resvecGP] = gmres(L, RHS, m, tol, maxcycles, diag(diag(L)));
relresvecGP = resvecGP ./ norm(RHS);
timeGP = toc();
disp(timeGP);
SOLGPRS = reshape(SOLGP, mx + 2, ny + 2);

% % BICGSTAB
% disp('BICG')
% tic
% [SOL4, flag4, res4 , ITER4, resvec4] = bicg(L, RHS, tol, maxit);
% relresvec4 = resvec4 ./ norm(RHS);
% time4 = toc();
% disp(time4)
% SOL4RS = reshape(SOL4, mx+2, ny+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
imagesc(flipud(uxy))
title('2D Poisson''s equation - analytical solution')
xlabel('x grid point')
ylabel('y grid point')
colorbar
set(gca, 'CLim', [0 100])
hold on;

figure(2)
% CHANGE PARAMETER 'SOL' HERE
% RS: rectangular grid
% SOLBNPRS: backslash, no preconditioning
% SOLBP: backslash, preconditioned system
% SOLPDRS: PD-GMRES, preconditioned system
% SOLLRS: LGMRES, preconditioned system
% SOLERS: GMRES-E, preconditioned system
imagesc(SOLLRS)
title('2D Poisson''s equation - numerical solution')
xlabel('x grid point')
ylabel('y grid point')
set(gca, 'YDir', 'Normal')
colorbar
hold on;

figure(3);
plot(mvec);
title('Evolution of restart parameter m')
hold on;

figure(4);
% GMRES(m) + prec.
semilogy(relresvecGP(1:m:size(relresvecGP, 1), 1), 'r', 'LineWidth', 2);
hold on;
% PD-GMRES + prec
semilogy(relresvecPD, 'm', 'LineWidth', 2);
hold on;
% LGMRES + prec
semilogy(relresvecL, 'g', 'LineWidth', 2);
hold on;
% GMRES-E + prec
semilogy(relresvecE, 'b', 'LineWidth', 2);
title('Relative residual norms');
xlabel('Number of restarts');
ylabel('||r_j||  / ||r_0 ||');
legend('GMRES(m)', 'PD-GMRES', 'LGMRES', 'GMRES-E');
hold on;