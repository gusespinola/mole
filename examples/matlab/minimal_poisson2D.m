% 2D Staggering example using a 2D Mimetic laplacian

clc
clear all
close all force

% You may want to change these 2 lines with 
% cd(FOLDER_TO_SCRIPT)
folder_to_cover = fullfile(getenv('MSCPATH'), ...
    'csrc-mole/examples/matlab/');
cd(folder_to_cover)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL PARAMETERS - DISCRETIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 4; % Order of accuracy
%h = 0.02; % Discretization size
%m = 1/h; % Vertical resolution
%n = 1/h; % Horizontal resolution
mx = 128;
ny = 128;
disp('k');
disp(k);
disp('mx');
disp(mx);
h = 1/mx;
hy = 1/ny;

% Robin parameters, alpha = 1 beta = 0 for Dirichlet;
alpha = 1; 
beta = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TURN THESE ON/OFF FOR SOME FEATURES (PLOT FIGURES, SAVE FILES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_preconditioned = 0;
plot_figures = 1;
save_experiment = 1;
save_linear_systems = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL PARAMETERS - PRECONDITIONING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment only once at a time
% preconditioner = 'none';
% preconditioner = 'Jacobi';
% preconditioner = 'Gauss_Seidel';
omega = 1.75; % Only useful for SOR
preconditioner = strcat('SOR_omega_', num2str(omega));
% preconditioner = 'ILU';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING ALREADY DISCRETIZED LINEAR SYSTEMS
% You may want to change these 2 lines with 
% cd(FOLDER_TO_LIN_SYS)
lin_sys_folder = fullfile(getenv('MSCPATH'), ...
    'csrc-mole/examples/matlab/linear_systems');
cd(lin_sys_folder)

linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
        '_m_', num2str(mx),'.mat')';

% if(strcmp(preconditioner, 'none')) == 1
%     % Will load existing non-preconditioned system
%     linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
%         '_m_', num2str(mx),'.mat')';
% else
%     % Will load existing preconditioned system
%     % First experiments will need to save this system before 
%     % its recycling
%     linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
%         '_m_', num2str(mx), '_', preconditioner, '.mat');
% end

load(linear_system);

cd('../');

CondNumber = condest(L);
disp('Condition Number before prec.');
disp(CondNumber);

fig7 = figure(7);
spy(L);
hold on;
figName = strcat('sparsity_k_', num2str(k), '_m_', num2str(mx),...
        '_none.fig');
savefig(fig7, figName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCRETIZATION - MUST CHANGE GLOBAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L2 = lap2D(k, mx, 1, ny, 1); % 2D Mimetic laplacian operator
% L = L2 + robinBC2D(k, mx, 1, ny, 1, alpha, beta);
% 
% RHS = zeros(mx+2, ny+2);
% 
% % Source term
% RHS(1, :) = 100; % Known value at the bottom boundary
% 
% RHS = reshape(RHS, [], 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYTICAL SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=0:h:1;
y=0:h:1;
[X, Y]=meshgrid(x, y);
uxy=zeros(size(X));
%tic
for n=1:100
    cn=(1-(-1)^n)/(sinh(n*pi)*n);
%     uxy=uxy+cn*sin(n*pi*X)...
%         .*((exp(n*pi*(Y))-exp(n*pi*(2.-Y))))...
%         /(1-exp(2*n*pi));
    uxy=uxy+cn*sin(n*pi*X).*sinh(n*pi*Y);
end
uxy=(200/pi).*uxy;
%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEFT-PRECONDITIONING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timePREC = 0;

disp(strcat('PRECONDICIONADOR = ', preconditioner));

% NO Preconditioning
disp(timePREC);

% Jacobi
if strcmp(preconditioner, 'Jacobi') == 1
    tic
    M_jacobi = (diag(diag(L)))^-1;
    L = M_jacobi * L;
    RHS = M_jacobi * RHS;
    timePREC = toc();
    disp(timePREC);
end

% Gauss-Seidel
if strcmp(preconditioner, 'Gauss_Seidel') == 1
    tic
    M_GS = (tril(L))^-1;
    L= M_GS*L;
    RHS = M_GS*RHS;
    timePREC = toc();
    disp(timePREC);
end

% SOR
if strcmp(preconditioner, strcat('SOR_omega_', num2str(omega))) == 1
    tic
    n=size(L);
    tStart=tic;
    omega=1;
    D=diag(diag(L));
    L=-(tril(L)-D);
    M_sor=(1/omega)*D-L;
    L=M_sor\L;
    RHS=M_sor\RHS;
    timePREC = toc();
    disp(timePREC);
end

% ILU
if strcmp(preconditioner, 'ILU') == 1
    tic
    setup.type = 'ilutp';
    setup.droptol = 1e-6;
    setup.udiag=1;
    
    [LT,UT] = ilu(sparse(L), setup);
    M=LT * UT;
    L=M\L;
    RHS=M\RHS;
    timePREC = toc();
    disp(timePREC);
end

% Estimating condition number of coefficient matrix AFTER preconditioning
CondNumber = condest(L);
disp('Condition Number after prec.');
disp(CondNumber);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRAL ANALYSYS - EIGENVALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Six largest-magnitude eigenvalues of coefficient matrix
% max_autovals = eigs(L);

% n smallest-magnitude eigenvalues of coefficient matrix
% nEigs = 10;
% min_autovals = eigs(L, nEigs, 'smallestabs');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NUMERICAL SOLUTION - BACKSLASH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('BACKSLASH');
tic;
SOLB = L\RHS;
timeB = toc();
SOLBRS = reshape(SOLB, mx+2, []);
disp(timeB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KRYLOV SUBSPACE METHODS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL PARAMETERS
tol = 1e-9;
maxit = 100;

% You may want to change these 2 lines with 
% cd(FOLDER_TO_SOLVERS)
solver_folder = fullfile(getenv('MSCPATH'), ...
    'krysbas-dev/test-prec-gmres/src/solvers');
cd(solver_folder)

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
kk = 3;
disp('LGMRES');
[SOLL, flagL, relresvecL, kdvecL, timeL] = ...
    lgmres(L, RHS, m, kk, tol, maxit);
disp(timeL);
SOLLRS = reshape(SOLL, mx + 2, ny + 2);

% GMRES-E
m = 27;
d = 3;
disp('GMRES-E');
[SOLE, flagE, relresvecE, kdvecE, timeE] = ...
    gmres_e(L, RHS, m, d, tol, maxit);
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
iternG = (iterG(1,1)-1)*m + iterG(1,2);
relresvecG = resvecG ./ norm(RHS);
relresvecGc = relresvecG(1:m:size(relresvecG, 1), 1);
timeG = toc();
disp(timeG);
SOLGRS = reshape(SOLG, mx + 2, ny + 2);

% GMRES(m) + PRECONDITIONER
% disp('PREC-GMRES(m), GAUSS-SEIDEL');
% tic;
% [SOLGP, flagGP, resGP, iterGP, resvecGP] = gmres(L, RHS, m, tol, maxcycles, tril(L));
%                                                 % Jacobi: M = diag(diag(L));
%                                                 % GS: M = tril(L);
% relresvecGP = resvecGP ./ norm(RHS);
% timeGP = toc();
% disp(timeGP);
% SOLGPRS = reshape(SOLGP, mx + 2, ny + 2);

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
if plot_figures == 1       
    figure(1)
    imagesc(flipud(uxy))
    title('2D Poisson''s equation - analytical solution')
    xlabel('x grid point')
    ylabel('y grid point')
    colorbar
    %set(gca, 'CLim', [0 100])
    set(gca, 'YDir', 'Normal')
    hold on;
    
    figure(2)
    % CHANGE PARAMETER 'SOL' HERE
    % RS: rectangular grid
    % SOLBRS: backslash
    % SOLPDRS: PD-GMRES
    % SOLLRS: LGMRES
    % SOLERS: GMRES-E
    imagesc(SOLLRS)
    title('2D Poisson''s equation - numerical solution')
    xlabel('x grid point')
    ylabel('y grid point')
    set(gca, 'YDir', 'Normal')
    colorbar
    hold on;
    
    fig3 = figure(3);
    plot(mvec);
    hold on;
    titleMVEC = strcat('m(j)\_k\_', num2str(k), '\_m\_', num2str(mx),...
        '\_prec\_', preconditioner);
    title(titleMVEC, 'FontSize', 16);
    figName = strcat('m(j)_k_', num2str(k), '_m_', num2str(mx),...
        '_prec_', preconditioner, '.fig');
    savefig(fig3, figName);
    
    fig4 = figure(4);
    % GMRES(m)
    semilogy(relresvecGc, 'r', 'LineWidth', 2);
    hold on;
    % PD-GMRES
    semilogy(relresvecPD, 'm', 'LineWidth', 2);
    hold on;
    % LGMRES
    semilogy(relresvecL, 'g', 'LineWidth', 2);
    hold on;
    % GMRES-E
    semilogy(relresvecE, 'b', 'LineWidth', 2);
    titleRRN = strcat('k = ', num2str(k), ', m = ', num2str(mx),...
        ', prec = ', preconditioner);
    title(titleRRN, 'FontSize', 16);
    xlim([0 100]);
    xticks([0 20 40 60 80 100]);
    xticklabels({'0','20','40','60','80','100'})
    ylim([1e-9 1]);
    yticks([1e-9 1e-6 1e-3 1]);
    yticklabels({'10^{-9}', '10^{-6}', '10^{-3}', '1'});
    xlabel('Number of restarts');
    ylabel('||r_j||  / ||r_0 ||');
    legendG = strcat('GMRES(30), cycles = ', num2str(iterG(1,1)), ...
        ', time = ', num2str(timeG), 's');
    legendL = strcat('LGMRES(27,3), cycles = ', num2str(size(relresvecL,1)), ...
        ', time = ', num2str(timeL), 's');
    legendE = strcat('GMRES-E(27,3), cycles = ', num2str(size(relresvecE,1)), ...
        ', time = ', num2str(timeE), 's');
    legendPD = strcat('PD-GMRES(30, \alpha_P=1, \alpha_D=1), cycles = ', num2str(size(mvec,1)), ...
        ', time = ', num2str(timePD), 's');
    legend(legendG, legendL, legendE, legendPD);
    pbaspect([4 3 1]); % or ([16 9 1])
    hold on;
    figName = strcat('relresvec_k', num2str(k), '_m_', num2str(mx),...
        '_prec_', preconditioner, '.fig');
    savefig(fig4, figName);
    
        figure(5)
        subplot(1,2,1)
        imagesc(uxy)
        title('2D Poisson''s equation - analytical solution')
        xlabel('x grid point')
        ylabel('y grid point')
        colorbar
        set(gca, 'Clim', [0 100])
        %set(gca, 'YDir', 'Normal')
        subplot(1,2,2)
        imagesc(SOLLRS)
        title('2D Poisson''s equation - numerical solution')
        xlabel('x grid point')
        ylabel('y grid point')
        colorbar
        %gca, 'Clim', [0 100])
        set(gca, 'YDir', 'Normal')
        hold on;
    
    if strcmp(preconditioner, "none") ~= 1
        figure(6)
        spy(L);
        hold on;
    end

    fig8 = figure(8);
    % GMRES(m)
    semilogy(relresvecGc, 'r', 'LineWidth', 2);
    hold on;
    titleRRN = strcat('k = ', num2str(k), ', m = ', num2str(mx),...
        ', prec = ', preconditioner);
    title(titleRRN, 'FontSize', 16);
    xlim([0 100]);
    xticks([0 20 40 60 80 100]);
    xticklabels({'0','20','40','60','80','100'})
    ylim([1e-9 1]);
    yticks([1e-9 1e-6 1e-3 1]);
    yticklabels({'10^{-9}', '10^{-6}', '10^{-3}', '1'});
    xlabel('Number of restarts');
    ylabel('||r_j||  / ||r_0 ||');
    legendG = strcat('GMRES(30), ciclos = ', num2str(iterG(1,1)), ...
        ', tiempo = ', num2str(timeG), 's');
    hold on;
    figName = strcat('relresvec_GMRESm_k', num2str(k), '_m_', num2str(mx),...
        '_prec_', preconditioner, '.fig');
    savefig(fig8, figName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE EXPERIMENT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_experiment == 1
    experiment_filename = strcat('experiment_k_',num2str(k),'_m_',...
        num2str(mx),'_',preconditioner,'_', datestr(now,30),'.mat');
    save(experiment_filename,...
        "k", "mx", "preconditioner", "CondNumber", ...
        "max_autovals", "min_autovals", ...
        "relresvecGc", "relresvecL", "relresvecE", "relresvecPD", ...
        "iterG", "iternG", "mvec",...
        "flagG", "flagE", "flagL", "flagPD", ...
        "timeB", "timeG", "timeL", "timeE", "timePD", "timePREC");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE NEW DISCRETIZED LINEAR SYSTEMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_linear_systems == 1 && strcmp(preconditioner, 'none') == 0
    new_linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
        '_m_', num2str(mx), '_', preconditioner, '.mat');
    save(new_linear_system, "L", "RHS");
end