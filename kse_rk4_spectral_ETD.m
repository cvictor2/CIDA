function kse_rk4_spectral_IF(N,T,lambda)
% Author: Adam Larios. Last modified on 2016-03-01.
% Solve the 1D Kuramoto-Sivashinsky equation on a periodic domain,
%   u_t + uu_x + lambda*u_xx + u_xxxx =0,    u(x,0) = u_0(x).
% Solution uses an explicit spectral method in space,
% RK4 time-stepping, 2/3's dealiasing, intergrating factor for linear term.

%Example:
% close all; kse_rk4_spectral_IF(512,4,0.001);

% mult_L = 2^(0);
% mult_u = 2^(0);
% mult_T = 2^(0);
% mult_lambda = 2^(0);

% mult_L = 2;
% mult_u = 1/8;
% mult_T = 16;
% mult_lambda = 4;

% mult_L = 1/(16*pi);
mult_L = 1;
mult_u = mult_L^(-3);
mult_T = mult_L^(4);
mult_lambda = mult_L^(-2)*6.748013066;  % Stable @ 6.748013065, Unstable @ 6.748013066

%% =========== Input Parameters ===========
close all;

% Plotting options
make_plots = 0; % 1= draw plot, else = don't draw plot
animate = 0; % 1 = yes, else = no
frames = 120; % number of frames in the animation.
params.nlTerm = 1;% 1 = use nonlinear term, else = no nonlinear term

left_endpoint = -16*pi*mult_L;
right_endpoint = 16*pi*mult_L;

%% === Pre-set options so you can just hit the "run" button in Matlab.
if ~(exist('N','var'))
    N = 2^10; % Number of gridpoints (uniform). Powers of 2 are most efficient for fft.
elseif mod(N,2)
    N = N+1;
    display('WARNING: Odd number of gridpoints not supported, adding one to make an even number.');
end
if ~(exist('T','var'))
%     T = 50*mult_T;
      T = 40;
end
if ~(exist('nu','var'))
    %lambda = 1*mult_lambda; 
    lambda = 1.64;
end

% N = 2^10, critical lambda = 1.64
% N = 2^11, critical lambda = 2.42
% N = 2^12, critical lambda = 3.25
% N = 2^13, critical lambda > 4.0


%% =========== Grid Set Up ===========
L = right_endpoint - left_endpoint;

% Physical space
dx = (right_endpoint - left_endpoint)/N; % Assume uniform grid.
x  = left_endpoint + (0:(N-1))*dx;

% Wave numbers k (we also multiply k by i for convenience).
params.ik   = 1i*[0:N/2-1 0 -N/2+1:-1]*(2*pi/L);
k_abs =      abs([0:N/2-1 0 -N/2+1:-1])*(2*pi/L);
k_sq =    [0:N/2-1 0 -N/2+1:-1].^2*(2*pi/L)^2;
k_qu =    [0:N/2-1 0 -N/2+1:-1].^4*(2*pi/L)^4;

params.dealias_modes = ceil(N/3):(N - floor(N/3)+1);

%% =========== Initial Conditions ===========
% Commment/Uncomment the initial condition you want

% % --- High and low frequency mixture.
% u = sin((2*pi/L)*x) + 0.25*cos(16*(2*pi/L)*x);
% % ---

% % --- Riemann Initial Conditions ---
% u = zeros(1,N); u(1:floor(2*N/3)) = 1;          %----_____
% u = zeros(1,N); u(floor(N/3):N) = 1;            %____-----
% u = zeros(1,N); u(floor(N/3):floor(2*N/3)) = 1;  %___---___
% % ---

% % --- Cubic interpolation from 0 to 1 ---
% ic_smooth_shock = @(x) (x<-1) + (0.25*x.^3-0.75*x+0.5).*((-1< x) - (1<x));
% u = ic_smooth_shock(x);
% % ---

% % --- Normalized Bump Function ---
% u = (exp(1./((x/mult_L).^2-1))).*((-1< (x/mult_L)) - (1<(x/mult_L)))/0.443993824382328;
% u(u == Inf) = 0; % 0*Inf = Inf, so zero-out those if they arise.
% u(isnan(u)) = 0; % 0*NaN = NaN, so zero-out those too.
% % ---

% % --- Square wave ---
% ic_square = @(x) ((-1< x) - (1<x));
% u = ic_square(x);
% % ---

% % --- Triangle wave ---
% ic_triangle = @(x) max(0,L/4-abs(x));
% u = ic_triangle(x);
% % ---

% % --- Calculus-type initial conditions ---
% u = -sin((2*pi/L)*x);
% u = max(0,1-abs(x)); % hat function
% u = x.^2;
% u = cos(x).*exp(sin(x));
% % ---

% % --- Random Fourier coefficients ---
% RNG_seed = randi(10000);
% % randn('state',RNG_seed); % Seed random number generator, to keep it from changing each time. (Matlab 2007 version)
% rng(RNG_seed);   % Seed random number generator, to keep it from changing each time. (Matlab 2014 version)
% display(sprintf('RNG seed value is %d',RNG_seed));
% Nmax = max(32,ceil(N/4)); % Nmax = number of non-zero Fourier coefficents
% z_real = randn(1,Nmax/2 - 1);
% z_imag = randn(1,Nmax/2 - 1);
% z = [0, z_real + 1i*z_imag, zeros(1,(N-Nmax)/2)]; 
% u = real(ifft([z conj(z(end:-1:1))])); % Impose anti-symmtery, and use ifft
% u = u/max(abs(u)); % Normalize maximum amplitude.
% % ---

% % --- Random noise ---
% u = 2*rand(1,N)-1;
% % ---

% % --- Good KSE initial condition ---
u = cos(x/16).*(1+sin(x/16));
% % ---

% u = u*mult_u;
u_hat = fft(u); % FFT of initial condition.

%% =========== Initialization and Setup ===========

% Compute a time-step using the advective CFL condition of initial data.
cfl_adv = 0.4;
% dt = cfl_adv*dx/(8.2*max(abs(u)));
dt = cfl_adv*dx/(50*max(abs(u)));

% % Reduce time step using viscous CFL if necessary.
% if (nu > 0)
%     cfl_visc = 0.25;
%     dt_visc = cfl_visc*dx^2/nu;
%     dt = min(dt,dt_visc);
% end

t = 0:dt:T;
nT_steps = length(t);

% Exponential pieces for integrating factor method
% E  = exp((lambda*k_sq-k_qu)*dt/2);
% E2 = exp((lambda*k_sq-k_qu)*dt);

inflation = 1.7*lambda^2; % This is roughly the amount of increase in the L^\infty norm
y_axis_max = 1.1*inflation*max(abs(u));
y_axis_min = -y_axis_max;


if make_plots == 1
    
    if animate == 1
        % Initialize a storage array for the solution.
        u_all = zeros(nT_steps,N);
        NN = N/2 + 1; % Only need first N/2+1 coefficients, due to symmetry.
        u_hat_all = zeros(nT_steps,NN);
        u_x_hat_all = zeros(nT_steps,NN);
    end
    
    % Set up for waterfall plot.
    wf_skip = floor(length(t)/60);
    t_wf = t(1:wf_skip:end);
    u_wf = zeros(length(t_wf),N);
    ti_wf = 1; % index
    
    energy = zeros(1,length(t));
    enstrophy = zeros(1,length(t));
    
    % Initialization for computing L^infty norms, etc.
    u_min = inf;
    u_max = -inf;
end

% Exponential pieces for integrating factor method
L = (lambda*k_sq-k_qu); % The linear operator on the right-hand side.
E  = exp(L*dt/2);
E2 = exp(L*dt);
% Eda  = exp(((1)*lambda*k_sq-k_qu)*dt/2);
% Eda2 = exp(((1)*lambda*k_sq-k_qu)*dt);

% Compute (exp(x) - 1)/x accurately.
% See Higham, "Accuracy and Stability of Numerical Algorithms", p. 20
% y = exp(x)
% if y == 1
%   f = 1;
% else
%   f = (y - 1)/log(y);
% end

% Compute dt*(exp(dt*L) - 1)/(dt*L) accurately for exponential time differencing.
exp_dtL = exp(dt*L);
for i = 1:N
    if (abs(L(i)) < eps(2)) % if zero
        ETD(i) = dt;
    else
        ETD(i) = dt*((exp_dtL(i) - 1)/log(exp_dtL(i)));
    end
end

u_min = 10^10;
u_max = -10^10;

Nhalf = N/2;

fh1 = figure;
set(fh1,'units','normalized','position',[0,0,1,1]); %[x y width height]
subplot(1,2,1);
p1 = loglog(abs(u_hat(1:Nhalf))/N);
axis([1,Nhalf,1e-25,1]);
title('FFT(u)');
get(p1,'ydata');

subplot(1,2,2);
p2 = plot(x,u);
axis([-16*pi,16*pi,y_axis_min,y_axis_max]);
get(p2,'ydata');

%% =========== Main Time-Stepping Loop ===========
counter_av = 0;

tic; % track wall-time
for ti = 1:nT_steps
    
    u_min = min([min(u),u_min]);
    u_max = max([max(u),u_max]);
    
%     if make_plots == 1
%         if animate == 1
%             u_all(ti,:) = u;
%             u_hat_all(ti,:) = u_hat(1:NN)/N;
%             u_x_hat_all(ti,:) = params.ik(1:NN).*u_hat(1:NN)/N;
%         end
%         if (mod(ti,wf_skip) == 1)
%             u_wf(ti_wf,:) = u;
%             ti_wf = ti_wf + 1;
%         end
%         
%         % % Compute L^2 norm in physical space using Riemann integration.
%         % energy(ti) = norm(u)*sqrt(dx);
%         % enstrophy(ti) = norm(u_x)*sqrt(dx);
%         
%         % Compute L^2 norm in spectral space using Parseval's law (more accurate).
%         energy(ti) = norm(u_hat)/N;
%         enstrophy(ti) = norm(params.ik.*u_hat)/N;
%     end
    
    u = real(ifft(u_hat));
    if (mod(ti,400)==0)
        title(sprintf('u(%1.3f)',t(ti)));
        set(p1,'ydata',abs(u_hat(1:Nhalf))/N);
        drawnow;
        
        title(sprintf('u(%1.3f)',t(ti)));
        set(p2,'ydata',u);
        drawnow;
    end
    
    
    % Euler with Exponential Time-Differencing (ETD)
    k1 = compute_rhs_hat(u_hat,params);
    u_hat = E2.*u_hat + ETD.*k1;
    
    % RK4 with integrating factor (IF):
%     k2 = compute_rhs_hat(E.*(u_hat  + 0.5*dt*k1),params);
%     k3 = compute_rhs_hat( E.*u_hat  + 0.5*dt*k2 ,params);
%     k4 = compute_rhs_hat(E2.*u_hat  +  dt*E.*k3 ,params);
%     u_hat = E2.*u_hat + (dt/6)*(E2.*k1 + 2*E.*(k2 + k3) + k4);
    
    if (t(ti) > 39)
       counter_av = counter_av + 1;
       u_hat_sum = u_hat_sum + u_hat/N;
    end


end
u_hat_sum/counter_av;
loglog(abs(u_hat(1:Nhalf))/N);

set(p1,'ydata',abs(u_hat(1:Nhalf))/N);

% figure
% loglog(abs(u_hat(1:Nhalf))/N);

% u_min
% u_max

wall_time = toc; % Stop recording wall time.
output_string = 'Done computing solution. Wall time for computation = %g. dt=%g. Time steps=%d. Step speed=%1.2e.';
display(sprintf(output_string,wall_time,dt,nT_steps,wall_time/nT_steps));

return;

%% =========== Plotting ============
if make_plots == 1
    
    m=3;n=2; % number of rows and columns for plots
    
    fh = figure;
    % Make full screen:
    set(fh,'units','normalized','position',[0,0,1,1]); %[x y width height]
    
    sc = 4; subplot(m,n,sc);
    plot(t,energy,'r','linewidth',2);
    legend('||u(t)||_{L^2}','location','best');
    axis([t(1),T,0,max(energy)])
    xlabel('t');
    title('Energy');
    
    sc = 6; subplot(m,n,sc);
    plot(t,enstrophy,'b','linewidth',2);
    legend('||u_x(t)||_{L^2}','location','best');
    axis([t(1),T,0,max(enstrophy)])
    xlabel('t');
    title('Enstrophy');
% plot(energy,enstrophy)
% axis([0, max(energy), 0, max(enstrophy)]);
    
    %colormap([0 0 0]);
    colormap copper;
    
    sc = 2; subplot(m,n,sc);
    waterfall(x,t_wf,u_wf); view(10,70);
    axis([left_endpoint, right_endpoint, 0, T, u_min, u_max]);
    xlabel x, ylabel t, zlabel u, grid off

    if animate == 1

        subplot(m,n,3);
        p2=scatter3(1:NN,real(u_hat_all(1,1:NN)),imag(u_hat(1,1:NN)),ones(1,NN),1:NN);
        axis([0,NN,min(min(real(u_hat_all))),max(max(real(u_hat_all))),min(min(imag(u_hat_all))),max(max(imag(u_hat_all)))]);
        title('Fourier coefficients of $u$','interpreter','latex');
        % view(0,0);
        xlabel('wave numbers');
        ylabel('re');
        zlabel('im');

        subplot(m,n,5);
        p3=scatter3(1:NN,real(u_x_hat_all(1,1:NN)),imag(u_hat(1,1:NN)),ones(1,NN),1:NN);
        axis([0,NN,min(min(real(u_x_hat_all))),max(max(real(u_x_hat_all))),min(min(imag(u_x_hat_all))),max(max(imag(u_x_hat_all)))]);
        title('Fourier coefficients of $u_x$','interpreter','latex');
        % view(90,0);
        xlabel('wave numbers');
        ylabel('re');
        zlabel('im');
        
        subplot(m,n,1);
        p1=plot(x,u_all(1,:),'b');
        axis([left_endpoint, right_endpoint, u_min, u_max]);
        title(sprintf('Solution u(%1.3f)',t(1)));
        xlabel('x');
        ylabel('u(x,t)');
        
        % Update data of plots to animate (faster than replotting).
        skip = 1+floor(nT_steps/frames);
        for ti = skip:skip:nT_steps
            drawnow;
            set(p1,'xdata',x','ydata',u_all(ti,:)');
            title(sprintf('Solution u(%1.3f)',t(ti)));
            set(p2,'ydata',real(u_hat_all(ti,:)),'zdata',imag(u_hat_all(ti,:)));
            set(p3,'ydata',real(u_x_hat_all(ti,1:NN)),'zdata',imag(u_x_hat_all(ti,1:NN)));
        end
    end
end


end % =========== End function burgers_rk4_spectral ============

function rhs_hat = compute_rhs_hat(u_hat,params)
% Compute the fft of the right-hand side of the equation: - uu_x

% Dealias by setting the middle Fourier coefficient to zero.
% This is where Matlab stores the highest frequencies.
u_hat(params.dealias_modes) = 0;

% Compute the product in physical space using dealiased versions.
u        = real(ifft(           u_hat));
u_x      = real(ifft(params.ik.*u_hat));
uu_x_hat = fft(u.*u_x); % Multiply, then go back to spectral space.

rhs_hat = -uu_x_hat;
% rhs_hat = (- uu_x_hat)./(1+(8/512)^2*k_sq);

end % =========== End function rhs_hat ============
