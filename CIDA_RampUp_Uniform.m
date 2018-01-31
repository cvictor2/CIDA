%function Chafee_Infante(N1,T,nui,alphai,graph)
%Chafee-Infante eq. u_t - nu*u_xx = u - alpha*u^3
%basic constants defined by problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solves Chafee-Infante equation using data assimilation on a uniform grid
%
%Meant to use loaded data from Data directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all;
close all;
clear all;
global alpha nu x N;
%% Load Variables
% load('Data\2018.01.23.01.24.03.575.mat')
if( exist('T','var') == 0)
    alpha = .1;
    nu = .1;
    dt = 0.01;
%     error = ;
    min_nodes = 4;
    L = 1;
    N = 2^12;
    mu = 5;
    seed = randi(10000);
    T = 200;
end

graph = true;
% graph = false;


%% Initialize other constants

dx = L/N;
x = linspace(0,L,N);
timesteps = ceil(T/dt);

%% Initial Conditions
% seed = randi(10000);
% fprintf('seed: %d\n', seed);
rng(seed);
a = (10^(0))*rand(1,N/4);
k = 1:1:N/4;
u_0 = zeros(1,N);
for i = 1:N/4
    u_0 = u_0 + (a(i)*sin(2*pi*k(i)*x));
end
u_0 = u_0/(norm(u_0)*sqrt(dx))*10^-3;


%% Plot Setup
%plot initial data and define plot
if(graph)
    figure
    subplot(1,2,1);
    h = plot(x,u_0);
    hold on;
%     h2 = plot(x,u_0);
    h2 = plot(x,zeros(1,N));
    y_amp = max([1.01/sqrt(alpha), 1e-2]);
    axis([0 L -y_amp y_amp]);
    set(h, 'XDataSource','x');
    set(h, 'YDataSource','u');
    set(h2, 'XDataSource','x');
    set(h2, 'YDataSource','v');
    title ( sprintf ('u(x ,%1.3f)',0 ));

%     drawnow;
    f = subplot(1,2,2);
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
    subplot(1,2,1)
    h3 = scatter(x_nodes,zeros(1,length(x_nodes)));
    set(h3, 'XDataSource','x_nodes');
    set(h3, 'YDataSource','zeros(1,length(x_nodes))');
end
v = zeros(size(u));
%             error_DA(1) = norm(u_0);
ui = u(2:N-1);
vi = v(2:N-1);
%% Main Loop
% close all
% f = figure;
% line([N/4 N/4], [0 inf]);
% for k = 1:timesteps
offset = 0;
while(max(abs(u))<.9/sqrt(alpha)&&t<T)
    %solves pde
    offset = offset +1;
    t = t+dt;
    umv = u - v;
    %error_DA(k) = sqrt(dx)*norm(umv);
%     if(max(error_DA(k)>1e50))
%         return
%     end
%     Ih_umv = interp1(x_nodes,umv(i_nodes),x)';
%     Ih_umv(1) = []; Ih_umv(end) = [];

    ui = B\(ui.*(ones(N-2,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
    
    u(2:N-1) = ui;
    
%% Graphing fft
if graph&&~ishghandle(f)&&~ishghandle(h)
    return
end
if (ishghandle(f)&&graph && mod(offset,100) == 0)
    subplot(1,2,2);
    u_hat = fft(u);
    loglog(abs(u_hat(1:N/2))/N);
    %         hold on;
    y1=get(gca,'ylim');
    line([N/4 N/4],[1e-20 1e10]);
    line([1e-0 N/4],[1e-15 1e-15]);
    axis([1,N/2,1e-20,1e5]);
    title ( sprintf ('Spectrum of u at time: %1.3f',t ));
    drawnow;
end
end
if t >= T
    close all;
    return
end
error_DA = zeros(1,timesteps-offset);

    %Error calculation
    %% Data Assimiliation
for( k = 1:timesteps-offset)
    
        
    t = t+dt;
    umv = u - v;
    error_DA(k) = sqrt(dx)*norm(umv);
    if(max(error_DA(k)>1e50))
        return
    end
    Ih_umv = interp1(x_nodes,umv(i_nodes),x)';
    Ih_umv(1) = []; Ih_umv(end) = [];

    ui = B\(ui.*(ones(N-2,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
    
    u(2:N-1) = ui;
    
    vi = B\(vi+dt*(alpha*(-2*vi-(vi.^3))+ mu*Ih_umv)); %Convex Splitting Method
    v(2:N-1) = vi;
    %     error_DA(k+1) = norm(u-v);
    
    
    %% Graphing
    %redraw plot
    if graph&&~ishghandle(f)&&~ishghandle(h)
        return
    end
%     if(max(u)> 0.9/sqrt(alpha)&&graph && mod(k,5) == 0
    if(graph && mod(k,20) == 0)
        subplot(1,2,1);
        refreshdata(h,'caller');
        refreshdata(h2,'caller');
        refreshdata(h3,'caller');
        refreshdata;
        title ( sprintf ('u(x ,%1.3f)',t ));
%         drawnow;
        %
        subplot(1,2,2);
        u_hat = fft(u);
        loglog(abs(u_hat(1:N/2))/N);
        %         hold on;
        y1=get(gca,'ylim');
        line([N/4 N/4],[1e-20 1e10]);
        line([1e-0 N/4],[1e-15 1e-15]);
        axis([1,N/2,1e-20,1e5]);
        title ( sprintf ('Spectrum of u at time: %1.3f',t ));
        drawnow;
    end
    
    
end
if(graph)
    figure
    semilogy(dt*(offset+1):dt:T, error_DA);
end

%% Save Data
%         date = datestr(now,'yyyy.mm.dd.HH.MM.SS.FFF');
%         save(['Data/' date '.mat'],'L','N','T','dt','alpha','nu','mu','seed','int_nodes','min_nodes','error');




