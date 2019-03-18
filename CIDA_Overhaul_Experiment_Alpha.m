function [alpha_nodes,all_min_nodes] = CIDA_Overhaul_Experiment_Alpha(mu_input,type)
%Chafee-Infante eq. u_t - nu*u_xx = u - alpha*u^3
%basic constants defined by problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solves Chafee-Infante equation using data assimilation on a uniform grid
%
%Meant to use loaded data from Data directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all;
% close all;
% clear all;
% param = struct; %All of the parameters for this problem: alpha, nu, min_nodes, L, N, mu

%% Load Variables
% load('Data\2018.01.23.01.24.03.575.mat')
% load('Data/Overhauled/2018.04.10.23.18.30_Initial_Ramp_Only_nu=1.0e-04_alpha=1_mu=100.mat')

graph = false;           %Display Graphs
save_check = true;     %Save data after ramp up and at end
% save_check = false;
date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');


%     alpha = 1;
    alpha = 150; % alpha is in [0.0001, 150]
    % Error in Fourier nodes is ~ 1e-12 for 0.0001.

     nu = 7.5e-6;
%     nu = 5e-5;
%     nu = 5e-4;
    %     nu is in [7.5e-6,0.01]
    
    L = 1;
    N = 2^12;
%     mu = 1000;
    mu = mu_input;
    min_nodes = max(ceil(.25*nu^-.5),5);
    
    seed = randi(10000);
%          seed = 3045
    T_ramp = 15;
    T = 50;
    %     param = struct('alpha',alpha,'nu',nu,'L',L,'mu',mu,'min_nodes',min_nodes,'seed',seed,'T',T);
    
    
    %% Initialize other constants
    
    dx = L/N;
    %x = linspace(0,L,N);
    x = 0 + dx*(0:(N-1));
    t_ramp = 0;
    
        dt = .001; %0.001;
%     dt = .001*10;
    timesteps = ceil(T/dt);
    
    
    %% Initial Conditions
    % seed = randi(10000);
    % fprintf('seed: %d\n', seed);
    u_0 = zeros(1,N);
   
    % dt = dx/(1/sqrt(alpha)); % CFL, Courant Fredrich Lewy number for advection

    
    
    u = u_0';
    ui = u(2:N);
    t=0;
    
    %% Initialize Param
    param = struct('alpha',alpha,'nu',nu,'L',L, 'N',N,'mu',mu,'min_nodes',min_nodes,'seed',seed,'T',T, 't', t, 'dt', dt, 'x', x, 'dx', dx, 'u_0', u_0, 't_ramp',t_ramp,'t_adj',t-t_ramp,'u_ramp',u_0, 'timesteps',timesteps);
    
    %% Ramp up
    offset = 0;
    u_solo = struct('u',u,'type',0);
    
    
%     if(save_check)
%         date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
%         save(sprintf('Data/Overhauled/%s_Initial_Ramp_Only_nu=%.1d_alpha=%.0d_mu=%.0d.mat',date_string,nu,alpha,mu),'param');
%     end
    vars = [];

%% Set up variables to track
offset = 0; %comment this out is you only want to run to time T, not time t_ramp + T
    %Type Descriptions:
    %   0 - Normal solution, no data assimilation
    %   1 - Uniform Grid Data Assimilation
    %   2 - Standard Car Data Assimilation
    %  -1 - Retro Car, old car method
    %   3 - Car with uniform grid along length (Hybrid Method)
    vars = [];
    vars = repelem(struct,1,1);
    i = 0;

        % Standard uniform grid setup
        if(type ==1)
            i = i+1;
            vars(i).type = 1;
            vars(i).v = zeros(size(u));
            vars(i).int_nodes = 10;
            vars(i).i_nodes = floor(linspace(1,N-1,vars(i).int_nodes));
            vars(i).x_nodes = x(vars(i).i_nodes);
            vars(i).error = zeros(1,timesteps-offset+1);
        else
            i = i+1;
            vars(i).type = 2;
            vars(i).v = zeros(size(u));
            vars(i).int_nodes = 30;
            vars(i).i_nodes = 1:vars(i).int_nodes;
            vars(i).x_nodes = x(vars(i).i_nodes);
            vars(i).velocity = 30*dx/dt;
            vars(i).direction = 1;
            vars(i).error = zeros(1,timesteps-offset+1);         
        end
    %Templates (Don't adjust these, just make copies if you want to use an
    %           instance)
%     % Standard uniform grid setup
%     i = i+1;
%     vars(i).type = 1;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes = min_nodes;
%     vars(i).i_nodes = floor(linspace(1,N-1,vars(i).int_nodes));
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).error = zeros(1,timesteps-offset+1);
%     
%     % Standard car setup
%     i = i+1;
%     vars(i).type = 2;
%     vars(i).v = zeros(size(u));
%     vars(i).int_nodes = 30;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).velocity = 10;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
% 
%     %Retro Car Configuration
%     i = i+1;
%     vars(i).type = -1;
%     vars(i).v = zeros(size(u));
%     vars(i).length_car = 10;
%     vars(i).int_nodes = vars(i).length_car;
%     vars(i).i_nodes = 1:vars(i).int_nodes;
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).velocity = vars(i).length_car;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
    
%     %Uniform Grid on Car Configuration
%     i = i+1;
%     vars(i).type = 3;
%     vars(i).v = zeros(size(u));
%     vars(i).length_car = 10; %Length of car (in gridpoints)
%                              %Set to N to cover entire grid
%     vars(i).int_nodes = 2;
%     vars(i).i_nodes = floor(linspace(1,vars(i).length_car, vars(i).int_nodes));
%     vars(i).x_nodes = x(vars(i).i_nodes);
%     vars(i).velocity = 10;
%     vars(i).direction = 1;
%     vars(i).error = zeros(1,timesteps-offset+1);
    clear i;

% trials = 5;
trials = 100;
retrials = 1;
% nu_nodes = logspace(log10(0.01),log10(7.5e-6),trials);
alpha_nodes = logspace(log10(150),log10(0.0001),trials);

% nu_nodes = [7.5e-6];
all_min_nodes = zeros(length(alpha_nodes),retrials);
seeds = all_min_nodes;
for(run = 1:trials)
    for(rerun = 1:retrials)
        seed = randi(10000);
        rng(seed);
        % Generate new initial data
        a = (10^(0))*randn(1,floor(N/5));
        k = 1:1:floor(N/5);
        u_0 = zeros(1,N);
        for i = 1:floor(N/5)
            u_0 = u_0 + (a(i)*sin(2*pi*k(i)*x/L));
        end
        u_0 = u_0/(norm(u_0)*sqrt(dx))*10^-3;
        
        seeds(run,rerun) = seed;
        t = 0;
        u = u_0';
        ui = u(2:N);
%         nu = nu_nodes(run);
        alpha = alpha_nodes(run);
        param.t = t;
        param.u = u;
%         param.nu = nu;
        param.alpha = alpha;
        for(i = 1:length(vars))
            vars(i).v = zeros(size(u));
        end
        if(vars(1).type==1)
            endpt1 = 2; endpt2 = 500;
        else
            endpt1 = 0; endpt2 = 500;
        end
        % end
        %% Ramp up
        A = (nu/(dx^2))*gallery('tridiag',N-1,1,-2,1) + (1+2*alpha)*speye(N-1);
        B = speye(N-1) - dt*A;
        
        
        while(max(abs(u))<.8/sqrt(alpha)&&t<T_ramp-dt)
            %solves pde
            offset = offset +1;
            t = t+dt;
            param.t = t;
%                     param.t_adj = t - t_ramp;
            
            %     umv = u - v;
            %error_DA(k) = sqrt(dx)*norm(umv);
            ui = B\(ui.*(ones(N-1,1)+dt.*(alpha.*(-2-(ui.^2))))); %Convex Splitting Method
            u(2:N) = ui;
            u_solo.u = u;
            %         u_f = abs(fft(u)/N);
            %     umv_error(offset) = sqrt(dx)*norm(u-v);
        end
        t_ramp = t;
        param.t_ramp = t_ramp;
        param.t_adj = t-t_ramp;
        u_ramp = u;
        param.u_ramp = u_ramp;
        
        
        %% Solve System
        % for j = 1:timesteps-offset
        % error = ones(length(vars),1);
        % error = 1;
        while(endpt2-endpt1>1)
            error = 1;
            result = floor((endpt1+endpt2)/2);
            if(endpt2 == 200 && endpt1 == 2)
                result = 3;
            end
            u = u_ramp;
            t = 0;
            for(i = 1:length(vars))
                vars(i).int_nodes = result;
                if(vars(i).type == 2)
                    vars(i).i_nodes = 1:vars(i).int_nodes;
                    vars(i).x_nodes = x(vars(i).i_nodes);
                    vars(i).x_nodes_o = vars(i).x_nodes;
                else
                    vars(i).i_nodes = floor(linspace(1,N-1,vars(i).int_nodes));
                    vars(i).x_nodes = x(vars(i).i_nodes);
                end

            end
            while(error>5e-14&&t<T)
                t = t+dt;
                param.t = t;
                param.t_adj = param.t_adj + dt;
                %     u - u_ramp;
                %     U = parallel.pool.Constant(u);
                %     U.Value;
                %update the v's
                for (k = 1:length(vars))
                    error = max(sqrt(dx)*norm(u-vars(k).v));
                    switch vars(k).type
                        case {2,-1, 3}
                            %update car nodes
                            vars(k) = Update_Overhaul(param, vars(k));
                    end
                    umv = u - vars(k).v;
                    Ih_umv = LinTerp(param, vars(k), umv);
                    vars(k).v(2:N) = B\(vars(k).v(2:N)+dt*(alpha*(-2*vars(k).v(2:N)-(vars(k).v(2:N).^3))+ (mu)*Ih_umv)); %Convex Splitting Method
                end
                %update u
                u(2:N) = B\(u(2:N).*(ones(N-1,1)+dt.*(alpha.*(-2-(u(2:N).^2))))); %Convex Splitting Method
                u_solo.u = u;
            end
            if(t>=T)
                endpt1 = result;
            else
                endpt2 = result;
            end
        end
        % timesErrors(trials) = t;
        all_min_nodes(run,rerun) = endpt2;
        
        % date = datestr(now,'yyyy.mm.dd.HH.MM.SS');
        if(save_check)
            newDate = datestr(now,'yyyy.mm.dd.HH.MM.SS');
            save(sprintf('Data/Overhauled/%s_Initial_Ramp_Var_Config_nu=%.1d_alpha=%.0d_mu=%.0d.mat',newDate,nu,alpha,mu),'param','vars');
            
            save(sprintf('Data/Overhauled/%s_Uniform_Trials_alpha(varies)_nu=%.0d_mu=%.0d.mat',date_string,nu,mu),'param','vars','all_min_nodes','alpha_nodes','seeds');
        end
    end
end
if(vars(1).type) == 1
    ttype = 'uniform';
else
    ttype = 'car';
end
save(sprintf('Results/%s_Uniform_alpha_Trials_type=%s_nu=%.0d_mu=%.0d_finished.mat',date_string,ttype,nu,mu),'param','vars','all_min_nodes','alpha_nodes','seeds');
end

function [lin] = LinTerp(p,v,umv)
if(length(v.i_nodes)>1)
    switch v.type
        case 1
            lin = interp1(v.x_nodes,umv(v.i_nodes),p.x)';
            lin(end) = [];
        case {2,-1, 3}
            
            [~,index] = max(abs(v.x_nodes(2:end)-v.x_nodes(1:end-1)));
            %         if(maximum*p.dx) > ceil((p.L/v.int_nodes)*p.dx)
            if(index~=1)
                group1 = struct('i_nodes',v.i_nodes(1:index),'x_nodes',v.x_nodes(1:index));
                lend1 = (group1.i_nodes(1)-2)*p.dx;
                rend1 = (group1.i_nodes(end)+2)*p.dx;
                group2 = struct('i_nodes',v.i_nodes(index+1:end),'x_nodes',v.x_nodes(index+1:end));
                lend2 = (group2.i_nodes(1)-2)*p.dx;
                rend2 = (group2.i_nodes(end)+2)*p.dx;
                
                if(group1.i_nodes(1)==1)
                    lin1 = interp1([group1.x_nodes,rend1],[umv(group1.i_nodes)',0],p.x,'linear',0)';
                else
                    lin1 = interp1([0,group1.x_nodes,rend1],[0,umv(group1.i_nodes)',0],p.x,'linear',0)';
                end
                if(group2.i_nodes(end)==p.N)
                    lin2 = interp1([lend2,group2.x_nodes],[0,umv(group2.i_nodes)'],p.x,'linear',0)';
                else
                    lin2 = interp1([lend2,group2.x_nodes,1],[0,umv(group2.i_nodes)',0],p.x,'linear',0)';
                end
                lin = lin1 + lin2;
                
                lin(end) = [];
                
            else
                %No gap in the car
                if(v.i_nodes(1) ==0)
                    lin = interp1(v.x_nodes(2:end),umv(v.i_nodes(2:end)),p.x,'linear',0)';
                else
                    lin = interp1(v.x_nodes,umv(v.i_nodes),p.x,'linear',0)';
                    
                    
                end
                lin(end) = [];
            end
        otherwise
            lin = zeros(p.N-1,1);
            
    end
else
    if(v.x_nodes~=0)
        switch v.type
            case 1
                lin = interp1([0,v.x_nodes,1], [0,umv(v.i_nodes)',0], p.x, 'linear',0)';
                lin(end) = [];
            case {-1,2}
                lin = zeros(p.N,1);
                lin(v.i_nodes) = umv(v.i_nodes);
                lin(end) = [];
            case 3
                lend = v.i_nodes - ceil(v.length_car/2);
                rend = v.i_nodes + ceil(v.length_car/2);
                if(lend<1)
                    m = unique([0, p.x(lend+p.N), 1,1+p.x(rend),2]);
                    linlin = interp1([m, 1+v.x_nodes],[zeros(length(m),1)',umv(v.i_nodes)'], [p.x,p.x+1],'linear',0)';
                    lin = linlin(1:p.N) + linlin(p.N+1:end);
                elseif(rend>p.N+1)
                    m = unique([0, p.x(lend), 1,1+p.x(rend-p.N),2]);
                    linlin = interp1([m, v.x_nodes],[zeros(length(m),1)',umv(v.i_nodes)'], [p.x,p.x+1],'linear',0)';
                    lin = linlin(1:p.N) + linlin(p.N+1:end);
                elseif(lend==1&&rend~=p.N+1)
                    lin = interp1([0,v.x_nodes, p.x(rend), 1],[0,umv(v.i_nodes),0,0], [p.x],'linear',0)';
                elseif(rend==p.N+1&&lend~=1)
                    lin = interp1([0,p.x(lend),v.x_nodes, 1],[0,0,umv(v.i_nodes),0], [p.x],'linear',0)';
                elseif(rend==p.N+1&&lend==1)
                    lin = interp1([0,v.x_nodes, 1],[0,umv(v.i_nodes),0],[p.x],'linear',0)';
                else
                    lin = interp1([0, p.x(lend),v.x_nodes, p.x(rend), 1],[0, 0,umv(v.i_nodes),0,0], [p.x],'linear',0)';
                    
                end
                lin(end) = [];
        end
%         lin(end) = [];
    else
        lin = 0.*umv;
        lin(end) = [];
    end
    % lin = 1;
end
end

function [var] = Update_Overhaul(p, v)
% var = v;
% v.x_nodes = unique(floor(mod(v.x_nodes + v.velocity*p.t_adj,1)./p.dx).*p.dx);
switch v.type
    case -1
        v.i_nodes = unique(mod(v.i_nodes+v.direction*v.velocity,p.N));
%         if(v.i_nodes(1)==0)
%             v.i_nodes(1) = [];
%             v.i_nodes = [v.i_nodes, p.N];
%         end
        v.x_nodes = p.x(v.i_nodes);
    case 2
%         v.x_nodes = unique(floor(mod(v.x_nodes + v.direction*v.velocity*p.dt,1)./p.dx).*p.dx);
        v.x_nodes = unique(floor(mod(v.x_nodes_o + v.direction*v.velocity*p.t_adj,1)./p.dx).*p.dx);

        [maximum,index] = max(v.x_nodes(2:end)-v.x_nodes(1:end-1));
        if(maximum > p.dx*v.int_nodes)
%             v.i_nodes(1:index) = v.i_nodes(1:index);
            v.i_nodes = (v.x_nodes+p.dx)./p.dx;

        elseif(maximum == p.dx*2) %Gap caused when crossing over boundary
            v.x_nodes(1:index) = v.x_nodes(1:index) - p.dx;
            v.i_nodes = (v.x_nodes+p.dx)./p.dx;
            
        else
            v.i_nodes = (v.x_nodes+p.dx)./p.dx;
        end
        v.i_nodes = sort(v.i_nodes);
        v.x_nodes = sort(v.x_nodes);
%         if(v.i_nodes(1)==0)
%             v.i_nodes(1) = [];
%             v.i_nodes = [v.i_nodes];
%             v.x_nodes(1) = [];
%             v.x_nodes = [v.x_nodes, 1];
%         end
    case 3
        v.x_nodes = unique(floor(mod(v.x_nodes + v.direction*v.velocity*p.dt,1)./p.dx).*p.dx);
        v.i_nodes = (v.x_nodes+p.dx)./p.dx;
        if(v.x_nodes>=p.N+1)
            v.x_nodes = [0];
            v.i_nodes = [1];
        end
%         if(v.i_nodes(1)==0)
%             v.i_nodes(1) = [];
%             v.i_nodes = [v.i_nodes, p.N];
%             v.x_nodes(1) = [];
%             v.x_nodes = [v.x_nodes, 1];
%         end
    otherwise
        %         var = v;
end
var = v;
end