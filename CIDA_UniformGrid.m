%function Chafee_Infante(N1,T,nui,alphai,graph)
%Chafee-Infante eq. u_t - nu*u_xx = u - alpha*u^3
%basic constants defined by problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solves Chafee-Infante equation using data assimilation on a uniform grid
%Meant to use loaded data from Data directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all;
close all;
global alpha nu x N;
% 
% N = 2^12;
% alpha = .01;
% alpha = 10^(-10);
% alpha = logalpha(z);
% nu = 0.000102;
% T = 100;
graph = true;
%     graph = false;

%% Initialize other constants
%     tol = 1e-5;
%     L = 1;
dx = L/N;
x = linspace(0,L,N);
%     dt = .01;

timesteps = ceil(T/dt);

%% Initial Conditions
% seed = randi(10000);
% fprintf('seed: %d\n', seed);
rng(seed);
u_0 = 0.01*randn(1,N);
u_0(1) = 0;
u_0(N) = 0;

%% Plot Setup
%plot initial data and define plot
if(graph)
    h = plot(x,u_0);
    hold on;
    % h2 = plot(x,u_0);
    h2 = plot(x,zeros(1,N));
    y_amp = max([1.01/sqrt(alpha), 1e-2]);
    axis([0 L -y_amp y_amp]);
    set(h, 'XDataSource','x');
    set(h, 'YDataSource','u');
    set(h2, 'XDataSource','x');
    set(h2, 'YDataSource','v');
end


% Convex Splitting Method
A = (nu/(dx^2))*gallery('tridiag',N-2,1,-2,1) + (1+2*alpha)*speye(N-2,N-2);
B = speye(N-2,N-2) - dt*A;

u = u_0';
t=0;
%% Data Assimilation
% mu = 50;
int_nodes = min_nodes;
i_nodes = floor(linspace(1,N,int_nodes));
x_nodes = x(i_nodes);
if(graph)
    h3 = scatter(x_nodes,zeros(1,length(x_nodes)));
    set(h3, 'XDataSource','x_nodes');
    set(h3, 'YDataSource','zeros(1,length(x_nodes))');
end
v = zeros(size(u));
error_DA = zeros(1,timesteps);
%             error_DA(1) = norm(u_0);
ui = u(2:N-1);
vi = v(2:N-1);
global pass;
pass = 0;
%% Main Loop
for k = 1:timesteps
    %solves pde
    
    t = t+dt;
    umv = u - v;
    error_DA(k) = sqrt(dx)*norm(umv);
    Ih_umv = [interp1(x_nodes,umv(i_nodes),x)]';
    Ih_umv(1) = []; Ih_umv(end) = [];

    ui = B\(ui.*(ones(N-2,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
    
    u(2:N-1) = ui;
    
    %Error calculation
    %% Data Assimiliation
    vi = B\(vi+dt*(alpha*(-2*vi-(vi.^3))+ mu*Ih_umv)); %Convex Splitting Method
    v(2:N-1) = vi;
    %     error_DA(k+1) = norm(u-v);
    
    
    %% Graphing
    %redraw plot
    if graph&&~ishghandle(h)
        return
    end
    if(graph && mod(k,10) == 0)
        refreshdata(h,'caller');
        refreshdata(h2,'caller');
        refreshdata(h3,'caller');
        title ( sprintf ('u(x ,%1.3f)',t ));
        drawnow;
        
    end
end
if(graph)
    figure
    semilogy(dt:dt:T, error_DA);
end

%% Save Data
%         date = datestr(now,'yyyy.mm.dd.HH.MM.SS.FFF');
%         save(['Data/' date '.mat'],'L','N','T','dt','alpha','nu','mu','seed','int_nodes','min_nodes','error');





