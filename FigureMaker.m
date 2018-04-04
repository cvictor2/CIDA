% close all;
% clear all;
% 
% load('Data/min_nodes_experiment_results_grid1.mat');
% alphas  = repelem(nu_nodes,4)';
% mins = reshape(all_min_nodes',[400,1]);
% % scatter(alphas,mins);
% % scatter(1:10,1:10)
% 
% %Data  
% % x = ones(1,12);
% % y = [1 1 2 3 3 3 4 5 6 2 2 2];
% %Engine
% [uxy, jnk, idx] = unique([alphas,mins],'rows');
% szscale = histc(idx,unique(idx));
% %Plot Scale of 25 and stars
% scatter(uxy(:,1),uxy(:,2),'sizedata',szscale*25)
% 
% hold on; 
% plot(uxy(:,1),uxy(:,1).^(-0.5)/4);
% 
% set(gca,'xscale','log')
% % set(gca,'yscale','log')
% 
% set(gca,'fontsize', 18);
% xlabel('Nu Values');
% ylabel('Minimum Nodes Required');
% title('Minimum Nodes vs. Nu Values');
% 
% load('Data/min_nodes_experiment_results_mu2.mat');
% mu_all = [1,5:5:200];
% trials = 2;
% x = [];
% y = [];
% unstablex = [];
% unstabley = [];
% for(i = 1:length(nu_nodes))
%     for(j = 1:trials)
%         for(k = 1:length(mu_all))
%             if(mu_data(((i-1)*trials + j),k)==1)
%                 x = [x,nu_nodes(i)];
%                 y = [y,mu_all(k)];
%             end
%             if(mu_data(((i-1)*trials + j),k)==-1)
%                 unstablex = [unstablex,nu_nodes(i)];
%                 unstabley = [unstabley,mu_all(k)];
%             end
%         end
%     end
% end
% figure
% [uxy, jnk, idx] = unique([x.',y.'],'rows');
% szscale = histc(idx,unique(idx));
% %Plot Scale of 25 and stars
% scatter(uxy(:,1),uxy(:,2),'filled','sizedata',szscale*10)
% % scatter(x,y);
% set(gca,'xscale','log')
% set(gca,'fontsize', 18);
% xlabel('Nu Values');
% ylabel('Mu Values');
% title('Convergent Mu Values for Sweeping Probe');
% M = nu_nodes.^(-1/2)/4;
% hold on;
% scatter(unstablex,unstabley,'r');
% plot(nu_nodes,M);
% trials = length(min_nodes
load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.03.03.14.47.51_MinNodesResults_Grid_mu=100_alpha=1_largescale.mat')
[nu_trials,trials] = size(all_min_nodes);
alphas  = repelem(nu_nodes,trials)';
mins = reshape(all_min_nodes',[nu_trials*trials,1]);
[uxy, jnk, idx] = unique([alphas,mins],'rows');
% szscale = histc(idx,unique(idx));
%Plot Scale of 25 and stars
scatter(nu_nodes,all_min_nodes(:,1))%,'sizedata',szscale*25)

hold on; 

load('C:\Users\colli\Desktop\MATLAB\UCARE\CIDA\Results\2018.02.23.15.45.53_MinNodesResults_car_mu=100_alpha=1_largescale.mat')
[nu_trials,trials] = size(all_min_nodes(:,1));
alphas  = repelem(nu_nodes,trials)';
% mins = reshape(all_min_nodes',[nu_trials*trials,1]);
mins = all_min_nodes;
[uxy, jnk, idx] = unique([alphas,mins],'rows');
% szscale = histc(idx,unique(idx));
scatter(nu_nodes,all_min_nodes(:,1),'.r')%,'sizedata',szscale*25)

plot(uxy(:,1),((uxy(:,1).^(-0.5)/4)));

set(gca,'xscale','log')
% set(gca,'fontsize',14);
% set(gca,'yscale','log')

set(gca,'fontsize', 18);
xlabel('\nu Values');
ylabel('Minimum Nodes Required (m_h)');
title(sprintf('Minimum Nodes For \\mu = %d',mu));
legend('Uniform Grid','Sweeping Probe','m_h = (1/4)\nu^{-1/2}');
