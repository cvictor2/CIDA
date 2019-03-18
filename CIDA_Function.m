function [param, vars] = CIDA_Function()
%function Chafee_Infante(N1,T,nui,alphai,graph)
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

graph = true;           %Display Graphs
save_check = false;     %Save data after ramp up and at end
reset_vars = true;      %Override loaded variable configurations with new definitions
    alpha = 1;
    nu = 7.5e-6;
    %     nu is in [7.5e-6,0.01]
    
    L = 1;
    N = 2^12;
    mu = 200;
    min_nodes = max(ceil(.25*nu^-.5),5);
    
    seed = randi(10000)
%          seed = 7263
         T_ramp = 10;
    T = 15;
    %     param = struct('alpha',alpha,'nu',nu,'L',L,'mu',mu,'min_nodes',min_nodes,'seed',seed,'T',T);
    
    
    %% Initialize other constants
    
    dx = L/N;
    %x = linspace(0,L,N);
    x = 0 + dx*(0:(N-1));
    t_ramp = 0;
    
    dt = .001; %0.001;
    timesteps = ceil(T/dt);
    
    
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

    
    
    u = u_0';
    ui = u(2:N);
    t=0;
    
    %% Initialize Param
    param = struct('alpha',alpha,'nu',nu,'L',L, 'N',N,'mu',mu,'min_nodes',min_nodes,'seed',seed,'T',T, 't', t, 'dt', dt, 'x', x, 'dx', dx, 'u_0', u_0, 't_ramp',t_ramp,'t_adj',t-t_ramp,'u_ramp',u_0, 'timesteps',timesteps);
    
    %% Ramp up
    offset = 0;
    u_solo = struct('u',u,'type',0);

    % Convex Splitting Method
    nOnes = ones(N-1, 1);
    C = -2*diag(nOnes, 0) + diag(nOnes(1:N-2), -1) + 1*diag(nOnes(1:N-2), 1);
%     C = sparse(C);
    A = (nu/(dx^2))*C + (1+2*alpha)*eye(N-1);
    B = eye(N-1) - dt*A;
    

    while(max(abs(u))<.8/sqrt(alpha)&&t<T_ramp-dt)
        
%         if (graph && mod(offset,200) == 0)
%             if(~Graphing_Overhaul(param,u_solo,[]))
%                 return
%             end
%             %     pause;
%         end
        
        %solves pde
        offset = offset +1;
        t = t+dt;
        param.t = t;
%         param.t_adj = t - t_ramp;

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
    
    if(save_check)
        date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
        save(sprintf('Data/Overhauled/%s_Initial_Ramp_Only_nu=%.1d_alpha=%.0d_mu=%.0d.mat',date_string,nu,alpha,mu),'param');
    end
    vars = [];


%% Set up variables to track
offset = 0; %comment this out is you only want to run to time T, not time t_ramp + T
% if(( exist('vars','var') == 0)||reset_vars)
    %Type Descriptions:
    %   0 - Normal solution, no data assimilation
    %   1 - Uniform Grid Data Assimilation
    %   2 - Standard Car Data Assimilation
    %  -1 - Retro Car, old car method
    %   3 - Car with uniform grid along length (Hybrid Method)
%     vars = [];
    vars = repelem(struct('type',0,'v',zeros(size(u)),'int_nodes',10,'i_nodes', [1:10],'x_nodes',x(1:10),'error',zeros(1,timesteps-offset+1),'velocity',0,'direction',0,'length_car',10),2,1);
    i = 0;
        % Standard uniform grid setup
    i = i+1;
    vars(i).type = 1;
    vars(i).v = zeros(size(u));
    vars(i).int_nodes = min_nodes;
    vars(i).i_nodes = floor(linspace(1,N-1,vars(i).int_nodes));
    vars(i).x_nodes = x(vars(i).i_nodes);
    vars(i).error = zeros(1,timesteps-offset+1);

    % Standard car setup
    i = i+1;
    vars(i).type = 2;
    vars(i).v = zeros(size(u));
    vars(i).int_nodes = 30;
    vars(i).i_nodes = 1:vars(i).int_nodes;
    vars(i).x_nodes = x(vars(i).i_nodes);
    vars(i).velocity = 10;
    vars(i).direction = 1;
    vars(i).error = zeros(1,timesteps-offset+1);

    
    
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
% end

% P = parallel.pool.Constant(param);
% P.Value.
%% Solve System
for j = 1:timesteps-offset
    t = t+dt;
    param.t = t;
    param.t_adj = param.t_adj + dt;
%     U = parallel.pool.Constant(u);
%     U.Value;
    %update the v's
    for (k = 1:length(vars))
        vars(k).error(j) = sqrt(dx)*norm(u-vars(k).v);
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
    
    show = 100;
% show = 1;
    if(j<=50)
        show = 1;
    end
    %     if(max(u)> 0.9/sqrt(alpha)&&graph && mod(k,5) == 0
%     if(graph && mod(j,show) == 0)
%         if(~Graphing_Overhaul(param,u_solo,vars))
%             return
%         end
%     end
end
%Get last error measurement
 for k = 1:length(vars)
        vars(k).error(end) = sqrt(dx)*norm(u-vars(k).v);
 end
 
 %% Plot Error
 figure
 legendInfo = cell(1,length(vars));
 for k = 1:length(vars)
     if(k==2)
         semilogy(0:dt:T,vars(k).error,'--','LineWidth',1);
     else
         semilogy(0:dt:T,vars(k).error,'LineWidth',1);
     end
     hold on;
     switch vars(k).type
         case 0
             text = 'u';
         case 1
             text = 'Uniform Grid';
         case 2 
             text = 'Sweeping Probe';
         otherwise
             text = 'Misc.';
     end
     legendInfo{k} = [text];
 end
legend(legendInfo);
xlabel('Time');
ylabel('L^2 Norm of u-v');
title('Error Plot for Different Methods');
% date = datestr(now,'yyyy.mm.dd.HH.MM.SS');
if(save_check)
    if(~exists(date_string))
        date_string = datestr(now,'yyyy.mm.dd.HH.MM.SS');
    end
    save(sprintf('Data/Overhauled/%s_Initial_Ramp_Var_Config_nu=%.1d_alpha=%.0d_mu=%.0d.mat',date_string,nu,alpha,mu),'param','vars');
end
end


function [lin] = LinTerp(p,v,umv)
switch v.type
    case 1
        lin = interp1(v.x_nodes,umv(v.i_nodes),p.x)';
        lin(end) = [];
    case {2,-1, 3}
        [maximum,index] = max(v.x_nodes(2:end)-v.x_nodes(1:end-1));
        if(maximum*p.dx) > ceil((p.L/v.int_nodes)*p.dx)
%             maximum
            group1 = struct('i_nodes',v.i_nodes(1:index),'x_nodes',v.x_nodes(1:index));
            group2 = struct('i_nodes',v.i_nodes(index+1:end),'x_nodes',v.x_nodes(index+1:end));
            if(index ==1)
                lin = interp1(group2.x_nodes,umv(group2.i_nodes),p.x,'linear',0)';
                lin(1) = umv(1);
            elseif(length(group2.i_nodes)==1)
                lin = interp1(group1.x_nodes,umv(group1.i_nodes),p.x,'linear',0)';
                lin(p.N-1) = umv(p.N-1);
            elseif(v.i_nodes(1) ==0) %0 is included
                lin1 = interp1(group1.x_nodes(2:end),umv(group1.i_nodes(2:end)),p.x,'linear',0)';
                lin2 = interp1(group2.x_nodes,umv(group2.i_nodes),p.x,'linear',0)';
                lin = lin1 + lin2;


            else %0 is not included
                lin1 = interp1(group1.x_nodes,umv(group1.i_nodes),p.x,'linear',0)';
                lin2 = interp1(group2.x_nodes,umv(group2.i_nodes),p.x,'linear',0)';
                lin = lin1 + lin2;
            end
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
% lin = 1;
end

function [var] = Update_Overhaul(p, v)
% var = v;
% v.x_nodes = unique(floor(mod(v.x_nodes + v.velocity*p.t_adj,1)./p.dx).*p.dx);
switch v.type
    case -1
        v.i_nodes = unique(mod(v.i_nodes+v.direction*v.velocity,p.N));
        if(v.i_nodes(1)==0)
            v.i_nodes(1) = [];
            v.i_nodes = [v.i_nodes, p.N];
        end
        v.x_nodes = p.x(v.i_nodes);
    case 2
        v.x_nodes = unique(floor(mod(v.x_nodes + v.direction*v.velocity*p.dt,1)./p.dx).*p.dx);
        [maximum,index] = max(v.x_nodes(2:end)-v.x_nodes(1:end-1));
        if(maximum > p.dx*v.int_nodes)
            v.i_nodes(1:index) = v.i_nodes(1:index);
        elseif(maximum == p.dx*2) %Gap caused when crossing over boundary
            v.x_nodes(1:index) = v.x_nodes(1:index) - p.dx;
            v.i_nodes = v.x_nodes./p.dx;
            
        else
            v.i_nodes = v.x_nodes./p.dx;
        end    
    case 3
        v.x_nodes = unique(floor(mod(v.x_nodes + v.direction*v.velocity*p.dt,1)./p.dx).*p.dx);
         v.i_nodes = v.x_nodes./p.dx;
         if(v.i_nodes(1)==0)
             v.i_nodes(1) = [];
             v.i_nodes = [v.i_nodes, p.N];
             v.x_nodes = p.x(v.i_nodes);
         end
    otherwise
        %         var = v;
end
var = v;
end
% 
% function [check] = Graphing_Overhaul(p,u,v)
% persistent grapher_var fig grapher_u saver
% if(p.t == 0)
%     grapher_var = [];
%     fig = [];
%     grapher_u = [];
%     saver = [];
% 
% end
% m = length(v);
% if(~isempty(fig)&&~isvalid(fig)&&isempty(grapher_u))
%     clear fig;
% end
% if(isempty(fig))%||~isvalid(fig))
% %     set(gca,'fontsize', 16);
% %     set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     fig = gcf;
% %     fig = figure('DefaultAxesFontSize',16);
% % Set up figure properties:
% % Enlarge figure to full screen.
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % Get rid of tool bar and pulldown menus that are along top of figure.
% % set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% %     set(gcf, 'Position', get(0, 'Screensize'));
% %         figure('units','normalized','outerposition',[0 0 1 1])
% 
% 
% 
% 
% else
%     if(~ishghandle(fig))
%         clear grapher_var fig grapher_u
% 
%         check  = false;
%         return
%     end 
% end
% if(gcf~=fig)
% %     set(findall(gcf,'-property','FontSize'),'FontSize',16)
% 
%     figure(fig)
%     
% end
% if(isempty(grapher_u))
%     grapher_u.graph = plot(p.x,u.u, 'LineWidth',1);
%     set(grapher_u.graph, 'YDataSource','u.u');
%     hold on;
% 
% else
%     refreshdata(grapher_u.graph,'caller');
% end
% if(isempty(grapher_var))
%     grapher_var = repelem(struct,m,1);
%     
%     for i = 1:m
%         switch v(i).type
%             case 0
%                 grapher_var(i).graph = plot(p.x,v(i).u);
%                 set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).u',i));
%                 grapher_var(i).DAgraph = [];
%                 hold on;
%             case 1
%                 if(i==1)
%                     grapher_var(i).graph = plot(p.x,v(i).v,'-g','LineWidth',1);
%                     set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).v',i));
%                     grapher_var(i).DAgraph = scatter(v(i).x_nodes,zeros(1,v(i).int_nodes),'ob');
%                     set(grapher_var(i).DAgraph, 'XDataSource',sprintf('v(%i).x_nodes',i));
%                     
%                     
%                 else
%                     grapher_var(i).graph = plot(p.x,v(i).v,'--r','LineWidth',1);
%                     set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).v',i));
%                     grapher_var(i).DAgraph = scatter(v(i).x_nodes,zeros(1,v(i).int_nodes),'+r');
%                     set(grapher_var(i).DAgraph, 'XDataSource',sprintf('v(%i).x_nodes',i));
%                     
%                 end
% 
%                
%                 
%                 hold on;
%             case {-1, 2, 3}
%                 grapher_var(i).graph = plot(p.x,v(i).v);
%                 set(grapher_var(i).graph, 'YDataSource',sprintf('v(%i).v',i));
%                 
%                 
%                 grapher_var(i).DAgraph = scatter(v(i).x_nodes,zeros(1,v(i).int_nodes),'+r');
%                 set(grapher_var(i).DAgraph, 'XDataSource',sprintf('v(%i).x_nodes',i));
%                 
%                 hold on;
%             otherwise
%                 
%         end
%     end
%     y_amp = max([1.01/sqrt(p.alpha), 1e-2]);
%     axis tight manual;
%     axis([0 p.L -y_amp y_amp]);
%     title ( sprintf ('u(x, %1.3f)',p.t ));
%     xlabel('X');
%     ylabel('U(X)');
%     if(m~=0)
%         legend('Reference Solution','Uniform Solution','Uniform Gridpoints','Probe Solution','Probe Gridpoints');
%     else
%          legend('Reference Solution');
%     end
% else
%     for i = 1:m
%         refreshdata(grapher_var(i).graph,'caller');
%         if(~isempty(grapher_var(i).DAgraph))
%             refreshdata(grapher_var(i).DAgraph,'caller');
%         end
%     end
%     title ( sprintf ('u(x, %1.3f)',p.t_adj ));
% 
%     
% end
% drawnow;
% % savecheck = true;
% % filename = 'carexample2.gif';
% % % Capture the plot as an image
% % frame = getframe(fig);
% % im = frame2im(frame);
% % [imind,cm] = rgb2ind(im,256);
% % 
% % % Write to the GIF File
% % if  isempty(saver)&&savecheck
% %     imwrite(imind,cm,filename,'gif','DelayTime',0.1,  'Loopcount',inf);
% %     saver = true;
% % elseif(savecheck)
% %     imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'WriteMode','append');
% % end
% check = true;
% return
% end