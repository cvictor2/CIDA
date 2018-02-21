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
global u v_g v_c v_o x_nodes_g x_nodes_c x_nodes_o u_f umv_error
global car grid optimal

%% Experiment Parameters
nu_trials = 100;
mu_trials = 41;
% mu_trials = 10;
% mu_nodes = 100:1:110;
trials = 3;
mu_nodes = [1, 5:5:200];%linspace(1,200,20)
nu_nodes = logspace(log10(0.01),log10(7.5e-6),100);

%Storage Variables
mu_data = zeros(nu_trials*trials,mu_trials);

% num_trials = 4;
% all_min_nodes = zeros(length(nu_nodes),num_trials);
for(trials_nu = 1:nu_trials)
    %% Load Variables
    % load('Data\2018.01.23.01.24.03.575.mat')
    %     if( exist('T','var') == 0)
    alpha = 1;
    %alpha = 0;
    %     nu = 7.5e-6;
    nu = nu_nodes(nu_trials);
    %     dt = 0.01;
    %     error = ;
%     min_nodes = 10;
    M = 0.25*nu^-0.5;
    min_nodes = max(ceil(M),5);
    L = 1;
    N = 2^12;
%     mu = 50;
    %      seed = 3552;
    T = 25;
    %     end
    
    graph = false;
    car = true;
    grid = true;
    optimal = false;
    opt_tol = 1e-5;
    % graph = false;
    
    % graph = false;
    
    tol = 5e-15;

    %% Initialize other constants
    
    dx = L/N;
    %x = linspace(0,L,N);
    x = 0 + dx*(0:(N-1));
    
    
    %% Initial Conditions
    for(trials_repeated = 1:trials)
        seed = randi(10000);
%         seed = 1263
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
        %% Data Assimilation
        global velocity direction pass
        if(car)
            int_nodes_c = min_nodes;
            i_nodes_c = 1:int_nodes_c;
            x_nodes_c = x(i_nodes_c);
            velocity = length(i_nodes_c);
            direction = 1;
            pass = 0;
        end
        if(grid)
            int_nodes_g = min_nodes;
            i_nodes_g = floor(linspace(1,N-1,int_nodes_g));
            x_nodes_g = x(i_nodes_g);
        end
        
        
        v_c = zeros(size(u));
        v_g = v_c;
        v_o = v_c;
        %             error_DA(1) = norm(u_0);
        ui = u(2:N);
        vi_g = v_g(2:N);
        vi_c = v_c(2:N);
        vi_o = v_o(2:N);
        
        u_f = abs(fft(u)/N);
        %% Ramp up
        offset = 0;
        while(max(abs(u))<.8/sqrt(alpha)&&t<T)
            offset = offset +1;
            t = t+dt;
            ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
            u(2:N) = ui;
        end
        if t >= T
            close all;
            return
        end
        
        for(trials_mu = 1:mu_trials)
            mu = mu_nodes(trials_mu);
            error_DA_c = ones(1,timesteps-offset);
            error_DA_g = ones(1,timesteps-offset);
            error_DA_o = ones(1,timesteps-offset);
%             error = 1;
            
            u_ramp = u;
            
%             ui = u_ramp(2:N);
            
            %             while(abs(endpt2 - endpt1)>1)
            
            error = 1;
            vi_o = zeros(1,N-1)';
            vi_g = zeros(1,N-1)';
            vi_c = zeros(1,N-1)';
            v_g = zeros(1,N)';
            v_c = zeros(1,N)';
            v_o = zeros(1,N)';

            ui = u_ramp(2:N);
            
%             int_nodes_g = floor((endpt2+endpt1)/2);
%             i_nodes_g = floor(linspace(1,N-1,int_nodes_g));
%             x_nodes_g = x(i_nodes_g);
            
%             opt_tol = (endpt2+endpt1)/2;
            %% Data Assimiliation
            k = 0;
            while(mu_data(trials_nu+trials_repeated-1,trials_mu) == 0 && k <=timesteps - offset)
                %             for k = 1:timesteps-offset
                k = k+1;
                
                t = t+dt;
                if(grid)
                    umv1 = u - v_g;
                    Ih_umv_g = interp1(x_nodes_g,umv1(i_nodes_g),x)';
                    Ih_umv_g(end) = [];
                    
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
                
                ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
                
                u(2:N) = ui;
                if(car)
                    vi_c = B\(vi_c+dt*(alpha*(-2*vi_c-(vi_c.^3))+ mu*Ih_umv_c)); %Convex Splitting Method
                    v_c(2:N) = vi_c;
                    error_DA_c(k+1) = sqrt(dx)*norm(u-v_c);
                end
                if(grid)
                    vi_g = B\(vi_g+dt*(alpha*(-2*vi_g-(vi_g.^3))+ mu*Ih_umv_g)); %Convex Splitting Method
                    v_g(2:N) = vi_g;
                    error_DA_g(k+1) = sqrt(dx)*norm(u-v_g);
%                     error_DA_g(k+1)
                    
                end
                if(optimal)
                    vi_o = B\(vi_o+dt*(alpha*(-2*vi_o-(vi_o.^3))+ mu*Ih_umv_o)); %Convex Splitting Method
                    v_o(2:N) = vi_o;
                    error_DA_o(k+1) = sqrt(dx)*norm(u-v_o);
                end
                if(error_DA_g(k+1)>error_DA_c(k+1)&&error_DA_g(k+1)>1e-14)
                    mu_data(trials_nu+trials_repeated-1,trials_mu) = 1;
                end
                if(error_DA_g(k+1) > 1e2||error_DA_c(k+1)>1e2)
                    mu_data(trials_nu+trials_repeated - 1,trials_mu) = -1;
                end
%                 if(error_DA_c(k+1)>1e2)
%                     mu_data(trials_nu+trials_repeated - 1,trials_mu) = -2;
%                 end
                
            end
%             close all;
%             plot(error_DA_g);
%             hold on;
%             plot(error_DA_c);
%             set(gca,'yscale','log')
%             drawnow;
%             mu
%             pause;
%             error = sqrt(dx)*norm(u-v_g);
%             if(error <tol)
%                 endpt2 = int_nodes_g;
%             else
%                 endpt1 = int_nodes_g;
%             end
            %             end
            %         end
            %         min_nodes = endpt2;
            % end
%             opt_tol = endpt2;
%             min_nodes = (int_nodes_g);
            %         all_min_nodes(nu_trials, trials_repeated) = min_nodes;
%             mu_data(trials_nu+trials_repeated - 1, trials_mu) = result;
        end
        date = datestr(now,'yyyy.mm.dd.HH.MM.SS.FFF');
        save(['Data/' date '.mat'],'L','N','T','dt','alpha','nu','mu','seed','int_nodes_g','min_nodes','error');
        
    end
%     date = datestr(now,'yyyy.mm.dd.HH.MM.SS.FFF');
    save('Data/min_nodes_experiment_results_mu1.mat','mu_data','nu_nodes');
end

% function [x_nodes, i_nodes] = updateNodes(prev_nodes)
% global x N direction velocity pass;
% % if(min(abs(v)) < 1)
% %     vel = 1;
% % else
% %     vel = 20;
% % end
% if(prev_nodes(end)+direction*velocity >= N-1)
%     direction = -1;
%     pass = 1;
%     i_nodes = [N-length(prev_nodes):N-1];
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
% end
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
