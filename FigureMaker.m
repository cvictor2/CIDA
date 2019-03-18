
mu = 1000;
% mu = 500;
% mu = 300;
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.07.14.07.22.03_Uniform_Trials_type=uniform_alpha=1_mu=300_finished.mat')
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.09.28.07.50.06_Uniform_Trials_type=uniform_alpha=1_mu=1000_finished.mat')
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.09.29.02.40.59_Uniform_Trials_type=uniform_alpha=1_mu=500_finished.mat')
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.09.29.12.18.36_Uniform_Trials_type=uniform_alpha=1_mu=50_finished.mat')
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Important\matlab.mat')
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Important\experiment_mu=1000_uniform.mat')
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Important\experiment_mu=500_uniform.mat')

min_nodes = max(all_min_nodes,[],1);
all_min_nodes = min_nodes';
nu_nodes = nu_nodes(1,:);
[nu_trials,trials] = size(all_min_nodes);
alphas  = repelem(nu_nodes,trials)';
mins = reshape(all_min_nodes',[nu_trials*trials,1]);
[uxy, jnk, idx] = unique([alphas,mins],'rows');
% szscale = histc(idx,unique(idx));
%Plot Scale of 25 and stars
h1 = scatter(nu_nodes,all_min_nodes(:,1))%,'sizedata',szscale*25)
hold on; 

% 
% x = nu_nodes;
% x_1 = nu_nodes(100);
% x_2 = nu_nodes(97); 
% y_1 = log(all_min_nodes(100)); 
% y_2 = log(all_min_nodes(97));
% y = exp((y_2 - y_1)/(x_2 - x_1)*(nu_nodes - x_2) + y_2);
% 
% x_1 = nu_nodes(98);
% x_2 = nu_nodes(61); 
% y_1 = log(all_min_nodes(98)); 
% y_2 = log(all_min_nodes(61));
% z = exp((y_2 - y_1)/(x_2 - x_1)*(nu_nodes - x_2) + y_2);
% plot(x,z);
% plot(x,y);
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.06.30.22.20.28_Uniform_Trials_nu(varies)_alpha=1_mu=500_finished.mat');
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.06.30.16.13.49_Uniform_Trials_nu(varies)_alpha=1_mu=300_finished.mat')

X = ones(length(nu_nodes),2);
X(:,2) = log(nu_nodes');
Y = log(all_min_nodes);
% X'XB = X'Y;
B = (X'*X)\(X'*Y);
b = B(1);
m = B(2);
st_dev = norm(Y - (b+m*X(:,2)))/sqrt(length(X)-1);
b
st_dev
% createFits;

% figure;
% plot(nu_nodes, all_min_nodes); hold on;
% set(gca,'yscale','log');
xx = linspace(min(log(nu_nodes)),max(log(nu_nodes)),1000);
h2 = loglog(exp(xx),exp(m*xx + b), '-','color',[0 0.5 0], 'LineWidth',1); hold on;
h3 = loglog(exp(xx),exp(m*xx + b+st_dev),'--','color',[0 0.5 0], 'LineWidth',1); hold on;
h4 = loglog(exp(xx),exp(m*xx + b-st_dev),'--','color',[0 0.5 0], 'LineWidth',1); hold on;
% axis([0 (nu_nodes(1)) 0 (120)]);


% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.07.03.03.22.17_Uniform_Trials_nu(varies)_alpha=1_mu=300_finished.mat'); %Velocity = 30
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.09.24.23.53.07_Uniform_Trials_type=car_alpha=1_mu=1000_finished.mat');
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.09.29.07.24.57_Uniform_Trials_type=car_alpha=1_mu=500_finished.mat')
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Important\experiment_mu=500_car.mat')
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Important\experiment_mu=1000_car.mat')

load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Important\experiment_mu=500_car_alt.mat')

min_nodes = max(all_min_nodes,[],1);
all_min_nodes = min_nodes';
nu_nodes = nu_nodes(1,:);

% set(gca,'xscale','log')
% set(gca,'fontsize',14);
% set(gca,'yscale','log')

set(gca,'fontsize', 18);
xlabel('\nu Values');
ylabel('Minimum Nodes Required (m_h)');
title(sprintf('Minimum Nodes For \\mu = %d',mu));
legend('Uniform Grid','Sweeping Probe','m_h = (1/4)\nu^{-1/2}');

% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.07.02.18.16.45_Uniform_Trials_nu(varies)_alpha=1_mu=300_finished.mat');
% load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Data\Overhauled\2018.05.26.00.53.45_Initial_Ramp_Var_Config_nu=7.5e-06_alpha=1_mu=100.mat')
[nu_trials,trials] = size(all_min_nodes(:,1));
alphas  = repelem(nu_nodes,trials)';
% mins = reshape(all_min_nodes',[nu_trials*trials,1]);
mins = all_min_nodes;
[uxy, jnk, idx] = unique([alphas,mins],'rows');
% szscale = histc(idx,unique(idx));
h5 = scatter(nu_nodes,all_min_nodes(:,1),'.r');%,'sizedata',szscale*25)
hold on;
% plot(uxy(:,1),((uxy(:,1).^(-0.5)/4)));

% X = ones(length(nu_nodes),2);
% X(:,2) = log(nu_nodes');
% Y = log(all_min_nodes);
% % X'XB = X'Y;
% B = (X'*X)\(X'*Y);
% b_1 = B(1);
% m_1 = B(2);
% st_dev = norm(Y - (b_1+m_1*X(:,2)))/sqrt(length(X)-1);
% % createFits;
% 
% % figure;
% % plot(nu_nodes, all_min_nodes); hold on;
% xx = linspace(min(log(nu_nodes)),max(log(nu_nodes)),1000);
% loglog(exp(xx),exp(m_1*xx + b_1)); hold on;
% loglog(exp(xx),exp(m_1*xx + b_1+st_dev)); hold on;
% loglog(exp(xx),exp(m_1*xx + b_1-st_dev)); hold on;
% axis([0 nu_nodes(1) 0 120]);
% 
set(gca,'xscale','log')
% % set(gca,'fontsize',14);
set(gca,'yscale','log')

xx = linspace(min(log(nu_nodes)),max(log(nu_nodes)),1000);
% h2 = loglog(exp(xx),exp((-1/2)*xx + log(1/2)), '-','color','cyan', 'LineWidth',1); hold on;


set(gca,'fontsize', 18);
xlabel('\nu Values');
ylabel('Minimum Nodes Required (m_h)');
title(sprintf('Minimum Nodes For \\mu = %d',mu));
txt = sprintf('m_h = %.4f\\nu^{%.4f}', exp(b),m);
txt2 = sprintf('m_h = %.4fe^{\\pm %.4f}\\nu^{%.4f}', exp(b),st_dev,m);

legend([h1, h2, h3, h5],{'Uniform Grid', txt, txt2, 'Sweeping Probe'});
axis([0 (nu_nodes(1)) 0 (120)]);

% legend('Uniform Grid','Sweeping Probe','Exponential Approx','Error Bounds (.95)');
