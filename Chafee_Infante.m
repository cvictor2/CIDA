%function Chafee_Infante(N1,T,nui,alphai,graph)
%Chafee-Infante eq. u_t - nu*u_xx = u - alpha*u^3
%basic constants defined by problem
% close all; clear all;
close all;
numAlpha = 10;
numTrials = 10;
tol = 1e-5;
logalpha = logspace(-10,0,numAlpha);
end_nodes = zeros(100,numTrials);
for(z = 1:numAlpha)
    % graph = true; %Check to display graphs
    global alpha nu x N r vel;
    
    N = 2^12;
    %     if mod(N,2)
    %         N = N+1;
    %         display('WARNING: Odd number of gridpoints not supported, adding one to make an even number.');
    %     end
    % alpha = .01;
    % alpha = 10^(-10);
    alpha = logalpha(z);
    nu = 0.000102;
    T = 50;
    %     graph = true;
    graph = false;
    
    %% Initialize other constants
    L = 1;
    dx = L/N;
    x = linspace(0,L,N);
    r = 1;
    vel = 10;
    % dt_visc = (dx^2)/(2*nu);
    % dt_adv = 0.5*dx/1; % 0.5*dx/(linear coefficient)
    % dt = min([dt_visc,dt_adv]);
    % dt = dt_adv;
    dt = .01;
    
    timesteps = ceil(T/dt);
    
    %% Initial Conditions
    % global u_exact;
    % exact(0);
    % u_0 = sin(3*x*2*pi) +.1*cos(15*x*2*pi);
    % u_0 = sin(x*pi);
    % u_0 = sin(.5*x) +.1*cos(15*x);
    % rng(1000);
    for(i = 1:numTrials)
        
        seed = randi(10000);
        % seed = 6701;
        % seed = 9735;
        %     seed = 949
        %     fprintf('seed: %d\n', seed);
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
        
        %Initialize Error for manufactured solution
        % error_MS = zeros(1,timesteps+1);
        %% Data Assimilation
        % mu = 1;
        mu = 50;
        % i_nodes = 1:10;
        % i_nodes = linspace(1,N,10);
        error = 1;
        min_nodes = 6;
        while(error>tol)
            min_nodes = min_nodes+1;
            int_nodes = min_nodes;
            i_nodes = floor(linspace(1,N,int_nodes));
            
            %     i_nodes =1:floor((N/int_nodes)):N;
            if i_nodes(end) ~= N
                i_nodes = [ i_nodes, N];
            end
            % i_nodes = [1,N];
            % i_nodes(end-1) = [];
            % i_nodes([2,3,6,7,8,11,12,14]) = [];
            
            
            % i_nodes(5:10) = [];
            % i_nodes = [1,1110,2616,3571,4096];
            % i_2 = floor((i_nodes(2:length(i_nodes)) + i_nodes(1:length(i_nodes)-1))./2);
            % i_nodes = unique([i_nodes,i_2]);
            
            % i_nodes = ind';
            % i_2 = floor((i_nodes(2:length(i_nodes)) + i_nodes(1:length(i_nodes)-1))./2);
            % i_nodes = unique([i_nodes,i_2]);
            
            x_nodes = x(i_nodes);
            if(graph)
                h3 = scatter(x_nodes,zeros(1,length(x_nodes)));
                set(h3, 'XDataSource','x_nodes');
                set(h3, 'YDataSource','zeros(1,length(x_nodes))');
            end
            v = zeros(size(u));
            error_DA = zeros(1,timesteps);
            % error_DA(1) = norm(u_0);
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
                %     if(mod(k,10)==0)%&&(abs(v(i_nodes(9)) - u(i_nodes(9)))) < tol)
                %         [x_nodes,i_nodes] = updateNodes(i_nodes,v(i_nodes));
                %         [length(i_nodes), error_DA(k)]
                %         scatter(x_nodes,zeros(1,length(x_nodes)));
                %     end
                %     Ih_umv = [zeros(1,i_nodes(1)-1),interp1(x_nodes,umv(i_nodes),x(i_nodes)),zeros(1,N-i_nodes(end))]';
                %     Ih_umv = [zeros(1,i_nodes(1)-1),interp1(x_nodes,umv(i_nodes),x(i_nodes)),zeros(1,N-i_nodes(end))]';
                Ih_umv = [interp1(x_nodes,umv(i_nodes),x)]';
                Ih_umv(1) = []; Ih_umv(end) = [];
                %     if t>20
                %         i = find(abs(u) < .1)
                % %         pause
                %     end
                %Manufactured Solution
                %     exact(t); %update exact solution u_exact
                %     ui = B\(ui+dt*(alpha*(-2*ui-(ui.^3))+ forcing())); %Convex Splitting Method with forcing term
                %     error_MS(k+1) = norm([0;u_exact;0] - u);
                
                
                %Solution Methods
                %     ui = B\(ui+dt*(alpha*(-2*ui-(ui.^3)))); %Convex Splitting Method
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
                %% Testing
                %     if(t>20&&t<20+dt)
                %     i = find(abs(u)<.1);
                %     ind = i;
                %     end
            end
            %     end
            if(graph)
                figure
                semilogy(dt:dt:T, error_DA);
            end
            error = error_DA(end);
        end
        end_nodes(z,i) = min_nodes;
        % fprintf('Number of nodes: %d\n\n',length(i_nodes));
        % semilogy(0:dt:T, error_MS);
        %end
        
        %% Save Data
        date = datestr(now,'yyyy.mm.dd.HH.MM.SS.FFF');
        save(['Data/' date '.mat'],'L','N','T','dt','alpha','nu','mu','seed','int_nodes','min_nodes','error');
    end
end
save('Data/results.mat','end_nodes','logalpha');





