%Set up vectors to hold the means, 
%Generate 3 graphs:
%   All of the velocities for 10:120
%   All of the velocites but cut out v<1
%   Log scaled means for each m_h
close all;
j = 0;
mn = zeros(12,1);
all_data = zeros(12,135);
leg = cell(12,1);
fig1 = figure;
fig2 = figure;
fig3 = figure;

%% 10
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.09.30.17.28.21_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=10_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 20
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.09.30.23.08.41_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=20_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 30
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.10.50.49_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=30_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 40
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.11.25.53_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=40_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 50
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.11.56.42_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=50_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 60
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.12.21.35_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=60_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 70
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.13.03.27_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=70_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 80
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.16.42.44_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=80_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 90
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.17.01.37_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=90_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 100
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.17.28.37_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=100_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 110
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.18.42.20_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=110_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);
%% 120
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.10.01.18.05.19_Velocity_Trials_nu(varies)_alpha=1_mu=300_int_nodes=120_finished.mat')
len = 0;
newTimes = [];
for(i=1:135)
    if(timesErrors(i)~=200&&i*param.dt/param.dx>=1)
        newTimes(len+1) = timesErrors(i);
    end
end
j = j+1;
all_data(j,:) = timesErrors';
% figure(fig1);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% figure(fig2);
% plot(velocities.*param.dt/param.dx,timesErrors, 'LineWidth',1); hold on; 
% axis([1 130 0 210]);
mn(j) = mean(newTimes);
leg{j} = sprintf('%d', vars(1).int_nodes);

%% Finish Graphs
figure(fig1)
% leg = cell;
% v = [1:12];
v = [2:2:12];
leg = {};
% leg = cell{};
for(i=1:length(v))
   plot(velocities(6:end).*param.dt/param.dx, all_data(v(i),(6:end)),'LineWidth',1); hold on; 
   leg{i} = sprintf('M = %d', v(i)*10);
end
set(gca, 'yscale','log');
set(gca, 'xscale','log');
legend(leg);
axis tight;
title('Convergence Rates for Probes of Different Lengths');
set(gca,'fontsize',18);
xlabel('Velocities ($v \cdot\frac{\Delta t}{\Delta x}$ )','Interpreter','latex');
ylabel('Time to Convergence');

figure(fig2);
v = [2:2:12];
leg = {};
% leg = cell{};
for(i=1:length(v))
   plot(velocities.*param.dt/param.dx, all_data(v(i),:),'LineWidth',1); hold on; 
   leg{i} = sprintf('M = %d', v(i)*10);
end
set(gca, 'yscale','log');
set(gca, 'xscale','log');
legend(leg);
axis tight;
title('Convergence Rates for Probes of Different Lengths');
set(gca,'fontsize',18);
xlabel('Velocities ($v \cdot\frac{\Delta t}{\Delta x}$ )','Interpreter','latex');
ylabel('Time to Convergence');
% ylabel('e^x', 'Interpreter','tex')



figure(fig3)
h1 = scatter([10:10:120],mn); hold on;
set(gca, 'yscale','log');
set(gca, 'xscale','log');
axis tight;
axis([8 130 2 40]);
% axis[(0 )]
xlabel('Number of Nodes');
ylabel('Time to Convergence');

X = ones(12,2);
X(:,2) = log([10:10:120]');
Y = log(mn);
% X'XB = X'Y;
B = (X'*X)\(X'*Y);
b = B(1);
m = B(2);
st_dev = norm(Y - (b+m*X(:,2)))/sqrt(length(X)-1);
xx = linspace(min(log([10:10:120])),max(log([10:10:120])),1000);
h2 = loglog(exp(xx),exp(m*xx + b), '-','color',[0 0.5 0], 'LineWidth',1); hold on;
h3 = loglog(exp(xx),exp(m*xx + b+st_dev),'--','color',[0 0.5 0], 'LineWidth',1); hold on;
h4 = loglog(exp(xx),exp(m*xx + b-st_dev),'--','color',[0 0.5 0], 'LineWidth',1); hold on;
txt = sprintf('T = %.3fM^{%.3f}', exp(b),m);
legend([h1, h2, h3],{'Mean Convergence Rates', txt, 'Std Dev'});
set(gca,'fontsize', 18);
title('Mean Convergence Rate for Probe Lengths');


