close all;
clear all;
% load('Data/comparison2.mat');
% alphas  = repelem(alpha_nodes,10)';
% mins = reshape(all_min_nodes',[200,1]);
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
% scatter(uxy(:,1),uxy(:,2),'filled','sizedata',szscale*25)
% set(gca,'xscale','log')
% set(gca,'fontsize', 18);
% xlabel('Alpha Values');
% ylabel('Minimum Nodes Required');
% title('Minimum Nodes vs. Alpha Values');

% load('Data/min_nodes_experiment_results3.mat');
% alphas  = repelem(nu_nodes,20)';
% mins = reshape(all_min_nodes',[200,1]);
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
% scatter(uxy(:,1),uxy(:,2),'filled','sizedata',szscale*25)
% set(gca,'xscale','log')
% set(gca,'fontsize', 18);
% xlabel('Alpha Values');
% ylabel('Minimum Nodes Required');
% title('Minimum Nodes vs. Alpha Values');

load('Data/min_nodes_experiment_results_grid1.mat');
alphas  = repelem(nu_nodes,4)';
mins = reshape(all_min_nodes',[400,1]);
% scatter(alphas,mins);
% scatter(1:10,1:10)

%Data  
% x = ones(1,12);
% y = [1 1 2 3 3 3 4 5 6 2 2 2];
%Engine
[uxy, jnk, idx] = unique([alphas,mins],'rows');
szscale = histc(idx,unique(idx));
%Plot Scale of 25 and stars
scatter(uxy(:,1),uxy(:,2),'sizedata',szscale*25)

hold on; 
plot(uxy(:,1),uxy(:,1).^(-0.5)/4);

set(gca,'xscale','log')
% set(gca,'yscale','log')

set(gca,'fontsize', 18);
xlabel('Nu Values');
ylabel('Minimum Nodes Required');
title('Minimum Nodes vs. Nu Values');

% hold on;
% load('Data/min_nodes_experiment_results10.mat');
% alphas  = repelem(nu_nodes,2)';
% mins = reshape(all_min_nodes',[200,1]);
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
% % figure
% scatter(uxy(:,1),uxy(:,2),'filled','sizedata',szscale*25)
% 
% 
% set(gca,'xscale','log')
% % set(gca,'yscale','log')
% set(gca,'fontsize', 18);
% xlabel('Nu Values');
% ylabel('Minimum Nodes Required');
% title('Minimum Nodes vs. Nu Values');
% figure;
% % load('Data/comparison.mat');
% u1 = u1./max(u1);
% u2 = u2./max(u2);
% plot(x,u1);
% hold on;
% plot(x,u2);
% legend(sprintf('alpha = %.3f',alpha1),sprintf('alpha = %.3f',alpha2));
% xlabel('X Values');
% ylabel('u(x,t)');
% title('U(x,t) for different values');
% t2 = 11.08;
% t1 = 14.78;
