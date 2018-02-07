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
global u v x_nodes u_f umv_error
%% Load Variables
% load('Data\2018.01.23.01.24.03.575.mat')
% if( exist('T','var') == 0)
alpha_nodes = logspace(-3,0,20);
all_min_nodes = zeros(20,10);
for alpha_trials = 1:20
    
    alpha = alpha_nodes(alpha_trials);
    %     alpha = .001;
    nu = .0001;
    %     dt = 0.01;
    %     error = ;
    min_nodes = 2;
    L = 1;
    N = 2^12;
    mu = 10;
    %     seed = randi(10000)
    %     seed = 6700;
    T = 100;
    % end
    %
%         graph = true;
    graph = false;
    
    
    %% Initialize other constants
    
    dx = L/N;
    %x = linspace(0,L,N);
    x = 0 + dx*(0:(N-1));
    
    
    %% Initial Conditions
    for trials = 1:10
        seed = randi(10000);
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
        %%Data Assimilation
        
        % mu = 50;
        % int_nodes = min_nodes;
        int_nodes = 10;
        i_nodes = floor(linspace(1,N-1,int_nodes));
        x_nodes = x(i_nodes);
        
        v = zeros(size(u));
        %             error_DA(1) = norm(u_0);
        ui = u(2:N);
        vi = v(2:N);
        u_f = abs(fft(u)/N);
        % umv_error = abs(u-v);
        % umv_error = ones(1,timesteps);
        %% Plot Setup
        %plot initial data and define plot
        if(graph)
            global h1 h2 h3 f1 e1
            figure
            subplot(1,2,1);
            h1 = plot(x,u_0);
            hold on;
            h2 = plot(x,zeros(1,N));
            h3 = scatter(x_nodes,zeros(1,length(x_nodes)));
            
            y_amp = max([1.01/sqrt(alpha), 1e-2]);
            axis([0 L -y_amp y_amp]);
            set(h1, 'XDataSource','x');
            set(h1, 'YDataSource','u');
            set(h2, 'XDataSource','x');
            set(h2, 'YDataSource','v');
            set(h3, 'XDataSource','x_nodes');
            %     set(h3, 'YDataSource','zeros(1,length(x_nodes))');
            title ( sprintf ('u(x ,%1.3f)',0 ));
            
            %     drawnow;
            %     f = subplot(1,2,2);
            subplot(1,2,2);
            f1 = loglog(u_f(1:N/2));
            hold on;
            y1=get(gca,'ylim');
            line([N/4 N/4],[1e-20 1e10],'Color','red');
            line([1e-0 N/4],[1e-15 1e-15],'Color','red');
            axis([1,N/2,1e-20,1e5]);
            set(f1, 'YDataSource','u_f(1:N/2)');
            title ( sprintf ('Spectrum of u at time: %1.3f',t ));
            
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
        while(max(abs(u))<.9/sqrt(alpha)&&t<T)
            %% Graphing fft
            if (graph && mod(offset,100) == 0)
                if(~Graphing(t))
                    return
                end
                %     pause;
            end
            
            %solves pde
            
            offset = offset +1;
            t = t+dt;
            umv = u - v;
            %error_DA(k) = sqrt(dx)*norm(umv);
            ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
            u(2:N) = ui;
            u_f = abs(fft(u)/N);
            %     umv_error(offset) = sqrt(dx)*norm(u-v);
        end
%         if t >= T
%             close all;
%             return
%         end
        error_DA = zeros(1,timesteps-offset);
        
        u_ramp = u;
        error = 1;
        int_nodes = 5;
        t_ramp = t;
        %     i_nodes = floor(linspace(1,N-1,int_nodes));
        %     x_nodes = x(i_nodes);
        
        %     v = zeros(size(u));
        %             error_DA(1) = norm(u_0);
        %     ui = u(2:N);
        %     vi = v(2:N);
        while(error>1e-10)
            t = t_ramp;
            int_nodes = int_nodes+1;
            i_nodes = floor(linspace(1,N-1,int_nodes));
            x_nodes = x(i_nodes);
            
            v = zeros(size(u));
            %             error_DA(1) = norm(u_0);
            ui = u(2:N);
            vi = v(2:N);
            %% Data Assimiliation
            j = 0;
            while(j<=timesteps - offset&& error>1e-10)
%             for j = 1:timesteps-offset
                j=j+1;
                t = t+dt;
                umv = u - v;
                %             error_DA(k) = sqrt(dx)*norm(umv);
                %             if(max(error_DA(k)>1e50))
                %                 return
                %             end
                Ih_umv = interp1(x_nodes,umv(i_nodes),x)';
                %     Ih_umv(1) = [];
                Ih_umv(end) = [];
                
                ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
                
                u(2:N) = ui;
                
                vi = B\(vi+dt*(alpha*(-2*vi-(vi.^3))+ mu*Ih_umv)); %Convex Splitting Method
                v(2:N) = vi;
                %     error_DA(k+1) = norm(u-v);
                %     umv_error(offset+k) = sqrt(dx)*norm(u-v);
                
                error = sqrt(dx)*norm(u-v);

                
                %% Graphing
                %redraw plot
                show = 20;
                if(j<=10)
                    show = 20;
                end
                %     if(max(u)> 0.9/sqrt(alpha)&&graph && mod(k,5) == 0
                if(graph && mod(j,show) == 0)
                    if(~Graphing(t))
                        return
                    end
                end
            end
        end
        min_nodes = int_nodes;
        all_min_nodes(alpha_trials,trials) = int_nodes;
        % if(graph)
        %     figure
        %     semilogy(dt*(offset+1):dt:T, error_DA);
        % end
        
        %% Save Data
        date = datestr(now,'yyyy.mm.dd.HH.MM.SS.FFF');
        save(['Data/' date '.mat'],'L','N','T','dt','alpha','nu','mu','seed','int_nodes','min_nodes','error');
    end
    save('Data/min_nodes_experiment_results.mat','all_min_nodes','alpha_nodes');
end
%     save('Data/min_nodes_experiment_results.mat','all_min_nodes','alpha_nodes');

function [x_nodes, i_nodes] = updateNodes(prev_nodes)
global x N direction velocity pass;
% if(min(abs(v)) < 1)
%     vel = 1;
% else
%     vel = 20;
% end
if(prev_nodes(end)+direction*velocity >= N-1)
    direction = -1;
    pass = 1;
    i_nodes = [N-length(prev_nodes):N-1];
    x_nodes = x(i_nodes);
elseif(prev_nodes(1)+direction*velocity <=1)
    direction = 1;
    i_nodes = [2:length(prev_nodes)+1];
    x_nodes = x(i_nodes);
else
    % n = min(find((abs(uv)==max(abs(uv)))));
    i_nodes = (prev_nodes+direction*velocity);
    x_nodes = x(i_nodes);
    % x_nodes = [x(prev_nodes),x(n)];
    % i_nodes = [prev_nodes,n];
    %     x_nodes = unique([x(prev_nodes),x(n)]);
    %     i_nodes = unique([prev_nodes,n]);
end
end
function [check] = Graphing(t)
global h1 h2 h3 f1 e1
global x u_f u v error_DA int_nodes N x_nodes umv_error
if(~ishghandle(h1)||~ishghandle(f1))%||~ishghandle(e1))
    check = false;
    return;
end
check = true;
subplot(1,2,1);
refreshdata(h1,'caller');
refreshdata(h2,'caller');
refreshdata(h3,'caller');
% refreshdata;
title ( sprintf ('u(x ,%1.3f)',t ));
%         drawnow;
%
subplot(1,2,2);
refreshdata(f1,'caller');
% refreshdata;
title ( sprintf ('Spectrum of u at time: %1.3f',t ));

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

