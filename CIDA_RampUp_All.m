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
global u v_g v_c v_o v_h1 x_nodes_g x_nodes_c x_nodes_o x_nodes_h1c x_nodes_h1g u_f umv_error 
global car grid optimal hybrid1
%% Load Variables
% load('Data\2018.01.23.01.24.03.575.mat')
if( exist('T','var') == 0)
    alpha = 1;
    %alpha = 0;
%     nu = 7.5e-6;
%     nu = 0.0001;
    nu = 1e-4;
%     nu = 6.25e-4;
%     nu = 1.5625e-04

    %     nu is in [7.5e-6,0.01]
%         dt = 0.02;
    %     error = ;
%     M = .25*nu^-.5;
%     M = 0.25*nu^-0.5
%     min_nodes = max(ceil(M),5);
%     min_nodes = 25;
    min_nodes = 15;
    L = 1;
    N = 2^12;
    mu = 100;
    car_mu = 1;
%    seed = randi(10000)
    seed = 1851
%      seed = 1263
    T = 50;
end

graph = true;
car = true;
grid = true;
optimal = true;
hybrid1 = false;
opt_tol = 1e-1;
% graph = false;

%% Initialize other constants

dx = L/N;
%x = linspace(0,L,N);
x = 0 + dx*(0:(N-1));
t_ramp = 0;

%% Initial Conditions
% seed = randi(10000);
% fprintf('seed: %d\n', seed);
rng(seed);
a = (10^(0))*randn(1,floor(N/5));
k = 1:1:floor(N/5);
u_0 = zeros(1,N);
for i = 1:floor(N/5)
    u_0 = u_0 + (a(i)*sin(2*pi*k(i)*x/L));
end
u_0 = u_0/(norm(u_0)*sqrt(dx))*10^-3;

% dt = dx/(1/sqrt(alpha)); % CFL, Courant Fredrich Lewy number for advection
dt = .01;
% dt <= dx/velocity.  Our max velocity is max value u ever takes on, namely
% 1/sqrt(alpha);
timesteps = ceil(T/dt);


% Convex Splitting Method
A = (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1) + (1+2*alpha)*speye(N-1);
B = speye(N-1) - dt*A;

u = u_0';
t=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Assimilation
global velocity direction pass
if(car)
%     int_nodes_c = floor(min_nodes/2);
%     int_nodes_c = min_nodes;
    int_nodes_c = 10;
    i_nodes_c = 1:int_nodes_c;
    x_nodes_c = x(i_nodes_c);
%     velocity = length(i_nodes_c);
%     velocity = 40;
    velocity =  20;
    direction = 1;
    pass = 0;
end
if(grid)
    int_nodes_g = min_nodes;
    i_nodes_g =[1         430         859        1201        1543        1795        2048        2302        2556        2890        3225        3660        4095];
    int_ndoes_g = length(i_nodes_g);
%     i_nodes_g = floor(linspace(1,N-1,int_nodes_g));
    x_nodes_g = x(i_nodes_g);
end
if(hybrid1)
    int_nodes_h1g = floor(min_nodes/2);
    i_nodes_h1g = floor(linspace(1,N-1,int_nodes_h1g));
    x_nodes_h1g = x(i_nodes_h1g);
    
    int_nodes_h1c = floor((min_nodes+1)/2);
    i_nodes_h1c = 1:int_nodes_h1c;
    x_nodes_h1c = x(i_nodes_h1c);
    velocity = length(i_nodes_h1c);
%     velocity =  20;
    direction = 1;
    pass = 0;
end
global ramp
ramp = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_c = zeros(size(u));
v_g = v_c;
v_o = v_c;
v_h1 = v_c;
%             error_DA(1) = norm(u_0);
ui = u(2:N);
vi_g = v_g(2:N);
vi_c = v_c(2:N);
vi_o = v_o(2:N);
vi_h1 = v_h1(2:N);
u_f = abs(fft(u)/N);
% umv_error = abs(u-v);
% umv_error = ones(1,timesteps);
%% Plot Setup
%plot initial data and define plot
if(graph)
    global h1 h2_g h2_c h2_o hybrid1_graph h3 h4 h5 f1 e1 hybrid1_grid
    figure
%     subplot(1,2,1);
    h1 = plot(x,u_0);
    hold on;
%     h2 = plot(x,zeros(1,N));
    if(grid)
        h2_g = plot(x,zeros(1,N));
        set(h2_g, 'XDataSource','x');
        set(h2_g, 'YDataSource','v_g');
        h3 = scatter(x_nodes_g,zeros(1,length(x_nodes_g)),'ob');
        set(h3, 'XDataSource','x_nodes_g');
%         set(h3, 'YDataSource','zeros(1,length(x_nodes_g))');
    end
    if(car)
        h2_c = plot(x,zeros(1,N));
        set(h2_c, 'XDataSource','x');
        set(h2_c, 'YDataSource','v_c');
        h4 = scatter(x_nodes_c,zeros(1,length(x_nodes_c)),'sr');
        set(h4, 'XDataSource','x_nodes_c');
%         set(h4, 'YDataSource','zeros(1,length(x_nodes_c))');
    end
    if(optimal)
        h2_o = plot(x,zeros(1,N));
        set(h2_o, 'XDataSource','x');
        set(h2_o, 'YDataSource','v_o');
    end
    if(hybrid1)
        hybrid1_graph =  plot(x,zeros(1,N));
        set( hybrid1_graph, 'XDataSource','x');
        set( hybrid1_graph, 'YDataSource','v_h1');
        hybrid1_grid = scatter([x_nodes_h1c;x_nodes_h1g],zeros(1,length(x_nodes_h1c)+length(x_nodes_h1g)),'sr');
        set(hybrid1_grid, 'XDataSource','[x_nodes_h1c;x_nodes_h1g]');
    end
    y_amp = max([1.01/sqrt(alpha), 1e-2]);
    axis([0 L -y_amp y_amp]);
    set(h1, 'XDataSource','x');
    set(h1, 'YDataSource','u');

    %     set(h3, 'YDataSource','zeros(1,length(x_nodes))');
    title ( sprintf ('u(x ,%1.3f)',0 ));
    
    %     drawnow;
    %     f = subplot(1,2,2);
%     subplot(1,2,2);
%     f1 = loglog(u_f(1:N/2));
%     hold on;
%     y1=get(gca,'ylim');
%     line([N/4 N/4],[1e-20 1e10],'Color','red');
%     line([1e-0 N/4],[1e-15 1e-15],'Color','red');
%     axis([1,N/2,1e-20,1e5]);
%     set(f1, 'YDataSource','u_f(1:N/2)');
%     title ( sprintf ('Spectrum of u at time: %1.3f',t-t_ramp ));
    
    %     subplot(2,2,[3 4])
    % %     e1 = plot(umv_error);
    %     e1 = semilogy(dt:dt:T, umv_error);
    %
    % %     ylim([1e-15 1e5])%     set(gca, 'YScale', 'log')
    %     set(e1, 'YDataSource', 'umv_error');
    %     title(sprintf('Data Assimilation Error at time: %1.3f',t));
    
    %
end


%% Main Loop
% close all
% f = figure;
% line([N/4 N/4], [0 inf]);
% for k = 1:timesteps
%% Ramp up
offset = 0;

while(max(abs(u))<.8/sqrt(alpha)&&t<T)
    %% Graphing fft
    if (graph && mod(offset,100) == 0)
        if(~Graphing(t,t_ramp))
            return
        end
        %     pause;
    end
    
    %solves pde
    
    offset = offset +1;
    t = t+dt;
%     umv = u - v;
    %error_DA(k) = sqrt(dx)*norm(umv);
    ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
    u(2:N) = ui;
    u_f = abs(fft(u)/N);
    %     umv_error(offset) = sqrt(dx)*norm(u-v);
end
offset = 0;
if t >= T
    close all;
    return
end
t_ramp = t;
% dt = 1.1*dt;
error_DA_c = zeros(1,timesteps-offset);
error_DA_g = zeros(1,timesteps-offset);
error_DA_o = zeros(1,timesteps-offset);
error_DA_h1 = error_DA_o;
error = 1;
ramp = true;
if(optimal)
%     i_nodes_o = find(abs(u)<1e-1);
%     i_nodes_o = uniquetol(i_nodes_o,10/max(i_nodes_o));
%     i_nodes_o = floor((i_nodes_o(1:end-1) + i_nodes_o(2:end))/2);
%     i_nodes_o = [i_nodes_o;1;N-1];
%     i_nodes_o = unique(i_nodes_o);
%     i_nodes_o(end) = [];
%     int_nodes_o = length(i_nodes_o);
% i_nodes_o = [1;N/2;N-1];
% i_nodes_o =i_nodes_g;
% i_nodes_o(3)=[]
% int_nodes_o = length(i_nodes_o);
i_nodes_o = [1:580,1050:N-1];
int_nodes_o = length(i_nodes_o);
x_nodes_o = x(i_nodes_o);
% int_nodes_o = length(i_nodes_o);
% x_nodes
%     x_nodes_0 = abs
    if(graph)
%         subplot(1,2,1);
        h5 = scatter(x_nodes_o,zeros(1,length(x_nodes_o)),'+m');
        set(h5, 'XDataSource','x_nodes_o');
        h2_o = plot(x,zeros(1,N));
        set(h2_o, 'XDataSource','x');
        set(h2_o, 'YDataSource','v_o');
    end
end

% pause;
%% Data Assimiliation
for k = 1:timesteps-offset
    
%     if(k==1000)
%         dt = 20*dt;
%         A = (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1) + (1+2*alpha)*speye(N-1);
%         B = speye(N-1) - dt*A;
%     end
    t = t+dt;
    if(grid)
        umv1 = u - v_g;
%         Grid Moves
%         if(k~=1)
%             i_nodes_g(2:end-1) = mod((i_nodes_g(2:end-1)+1),N-1);
%             i_nodes_g = unique(i_nodes_g);
%             if(i_nodes_g(1) == 0)
%                 i_nodes_g(1:2) = [1,2];
%             end
%             x_nodes_g = x(i_nodes_g);
%         end
        Ih_umv_g = interp1(x_nodes_g,umv1(i_nodes_g),x)';
        Ih_umv_g(end) = [];
        
    end
    if(hybrid1)
        umvh1 = u - v_h1;
        Ih_umv_h1g = interp1(x_nodes_h1g,umvh1(i_nodes_h1g),x)';
        Ih_umv_h1g(end) = [];   
        if(k~=1)
            [x_nodes_h1c, i_nodes_h1c] = updateNodes(i_nodes_h1c);
        end
%         umv2 = u - v_c;
        %         Ih_umv = interpl(x_nodes_c,umv2(i_nodes_c),
        Ih_umv_h1c = [zeros(1,i_nodes_h1c(1)-1),interp1(x_nodes_h1c,umvh1(i_nodes_h1c),x(i_nodes_h1c)),zeros(1,N-i_nodes_h1c(end))]';
        Ih_umv_h1c(1) = [];
    end
    if(optimal)
        umv3 = u - v_o;
        if(mod(k,50)==0&&max(abs(umv3))>opt_tol)
%             if(max(abs(umv3))>1e-5)
%                 [~,newnode] = max(abs(umv3));
%                 i_nodes_o = [i_nodes_o;newnode];
%                 i_nodes_o = unique(i_nodes_o);
%                 int_nodes_o = length(i_nodes_o);
%                 x_nodes_o = x(i_nodes_o);
%                 if(graph)
% %                     subplot(1,2,1);
%                     h5 = scatter(x_nodes_o,zeros(1,length(x_nodes_o)),'+m');
%                     set(h5, 'XDataSource','x_nodes_o');
%                 end
%             end
                
        end
        Ih_umv_o = interp1(x_nodes_o,umv3(i_nodes_o),x)';
        Ih_umv_o(end) = [];
        
    end
    if(car)
        if(k~=1)
            [x_nodes_c, i_nodes_c] = updateNodes(i_nodes_c);
        end
        umv2 = u - v_c;
        %         Ih_umv = interpl(x_nodes_c,umv2(i_nodes_c),
        Ih_umv_c = [zeros(1,i_nodes_c(1)-1),interp1(x_nodes_c,umv2(i_nodes_c),x(i_nodes_c)),zeros(1,N-i_nodes_c(end))]';
        Ih_umv_c(1) = [];

    end
    %     umv = u - v;
    %     error_DA_c(k) = sqrt(dx)*norm(umv1);
    %     error_DA_c(k) = sqrt(dx)*norm(umv2);
    
    %     if(max(error_DA(k)>1e50))
    %         return
    %     end
    %     Ih_umv = interp1(x_nodes,umv(i_nodes),x)';
    %     Ih_umv(1) = [];
    % Ih_umv(end) = [];
    
    ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
    
    u(2:N) = ui;
    if(car)
        vi_c = B\(vi_c+dt*(alpha*(-2*vi_c-(vi_c.^3))+ (car_mu*mu)*Ih_umv_c)); %Convex Splitting Method
        v_c(2:N) = vi_c;
        error_DA_c(k+1) = sqrt(dx)*norm(u-v_c);
    end
    if(grid)
        vi_g = B\(vi_g+dt*(alpha*(-2*vi_g-(vi_g.^3))+ mu*Ih_umv_g)); %Convex Splitting Method
        v_g(2:N) = vi_g;
        error_DA_g(k+1) = sqrt(dx)*norm(u-v_g);
    end
    if(hybrid1)
        vi_h1 = B\(vi_h1+dt*(alpha*(-2*vi_h1-(vi_h1.^3))+ (mu)*(car_mu*Ih_umv_h1c + Ih_umv_h1g))); %Convex Splitting Method
        v_h1(2:N) = vi_h1;
        error_DA_h1(k+1) = sqrt(dx)*norm(u-v_h1);
    end
    if(optimal)
        vi_o = B\(vi_o+dt*(alpha*(-2*vi_o-(vi_o.^3))+ mu*Ih_umv_o)); %Convex Splitting Method
        v_o(2:N) = vi_o;
        error_DA_o(k+1) = sqrt(dx)*norm(u-v_o);
    end
    %     vi_g = B\(vi_g+dt*(alpha*(-2*vi_g-(vi_g.^3))+ mu*Ih_umv)); %Convex Splitting Method
    %     v_g(2:N) = vi_g;
    %     error_DA(k+1) = norm(u-v);
    %     umv_error(offset+k) = sqrt(dx)*norm(u-v);
    
    
    %% Graphing
    %redraw plot
    show = 20;
    if(k<=10)
        show = 20;
    end
    %     if(max(u)> 0.9/sqrt(alpha)&&graph && mod(k,5) == 0
    if(graph && mod(k,show) == 0)
        if(~Graphing(t,t_ramp))
            return
        end
    end
    
    
end
if(graph)
    figure
    if(grid)
%         subplot(1,3,1)
        semilogy(dt*(offset):dt:T, error_DA_g);
    end
    hold on;
    if(car)
%         subplot(1,3,2)
        semilogy(dt*(offset):dt:T, error_DA_c);
    end
    hold on;
    if(optimal)
%         subplot(1,3,3)
        semilogy(dt*(offset):dt:T, error_DA_o);
%         int_nodes_o
    end
end

%% Save Data
%         date = datestr(now,'yyyy.mm.dd.HH.MM.SS.FFF');
%         save(['Data/' date '.mat'],'L','N','T','dt','alpha','nu','mu','seed','int_nodes','min_nodes','error');

function [x_nodes, i_nodes] = updateNodes(prev_nodes)
global x N direction velocity pass;
if(pass ==1)
    i_nodes = 1:length(prev_nodes);
    x_nodes = x(i_nodes);
    pass = 0;
elseif(prev_nodes(end)+direction*velocity >= N-1)
    pass = 1;
    i_nodes = N-length(prev_nodes):N-1;
    x_nodes = x(i_nodes);
else
    i_nodes = prev_nodes+direction*velocity;
    x_nodes = x(i_nodes);
    
    
        
end
% if(min(abs(v)) < 1)
%     vel = 1;
% else
%     vel = 20;
% end
% if(prev_nodes(end)+direction*velocity >= N-1)
% %     direction = -1;
% %     pass = 1;
% %     i_nodes = [N-length(prev_nodes):N-1];
% %     x_nodes = x(i_nodes);
% %     direction = 1;
%     pass = 1;
%     i_nodes = 1:length(prev_nodes);
%     x_nodes = x(i_nodes);
% elseif(prev_nodes(1)+direction*velocity <=1)
%     direction = 1;
%     i_nodes = [2:length(prev_nodes)+1];
%     x_nodes = x(i_nodes);
% else
%     % n = min(find((abs(uv)==max(abs(uv)))));
%     i_nodes = (prev_nodes+direction*velocity);
%     x_nodes = x(i_nodes);
%     % x_nodes = [x(prev_nodes),x(n)];
%     % i_nodes = [prev_nodes,n];
%     %     x_nodes = unique([x(prev_nodes),x(n)]);
%     %     i_nodes = unique([prev_nodes,n]);
% end
end
function [check] = Graphing(t,t_ramp)
global h1 h2_g h2_c h2_o hybrid1_graph h3 h4 h5 f1 e1 hybrid1_grid
global grid car optimal ramp hybrid1
global x u_f u v error_DA int_nodes N x_nodes_c x_nodes_g x_nodes_o x_nodes_h1c x_nodes_h1g umv_error v_g v_c v_o v_h1
if(~ishghandle(h1))%||~ishghandle(f1))%||~ishghandle(e1))
    check = false;
    return;
end
check = true;
% subplot(1,2,1);
refreshdata(h1,'caller');
if(grid)
    refreshdata(h2_g,'caller');
    refreshdata(h3,'caller');
end
if(car)
    refreshdata(h2_c,'caller');
    refreshdata(h4,'caller');
end
if(optimal && ramp)
    refreshdata(h5,'caller');
    refreshdata(h2_o,'caller');
end
if(hybrid1)
    refreshdata(hybrid1_graph,'caller');
    refreshdata(hybrid1_grid,'caller');
end
% refreshdata;
title ( sprintf ('u(x ,%1.3f)',t-t_ramp ));
%         drawnow;
%
% subplot(1,2,2);
% refreshdata(f1,'caller');
% % refreshdata;
% title ( sprintf ('Spectrum of u at time: %1.3f',t-t_ramp ));

% subplot(2,2,[3 4]);
% refreshdata(e1,'caller');
% refreshdata;
% title ( sprintf ('Data Assimilation Error at time: %1.3f',t ));
%         u_hat = fft(u);
%         u_hat = fft(u);
%         loglog(abs(u_hat(1:N/2))/N);
%         hold on;
%         y1=get(gca,'ylim');
% %         line([N/4 N/4],[1e-20 1e10],'Color','red');
%         line([1e-0 N/4],[1e-15 1e-15],'Color','red');
%         axis([1,N/2,1e-20,1e5]);
%         title ( sprintf ('Spectrum of u at time: %1.3f',t ));
drawnow;
end



