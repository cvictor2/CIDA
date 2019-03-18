%% Create dummy figures for demonstration and save them in current folder 
% The first figure does not contain a legend
% f = figure; 
% surf(peaks) 
% hgsave(f,'one.fig') 
% close(f) 
% % The second figure contains a legeng and requires slightly different
% % treatment
% f = figure; 
% t = 0:pi/50:10*pi; 
% plot(t, sin(t));
% legend(gca, 'sinx')
% hgsave(f,'two.fig') 
% close(f) 
%% New figure 
f = figure; 
a1 = subplot(2,2,1); 
a2 = subplot(2,2,2); 
a3 = subplot(2,2,3); 
a4 = subplot(2,2,4); 

%% Copy first existing figure on first subplot 
% Open exiting figure 
f_c = openfig('ExampleGraphs/Part1.fig'); 
% Identify axes to be copied 
axes_to_be_copied = findobj(f_c,'type','axes'); 
% title = findobj(f_c,'type','title');
% Identify the children of this axes 
chilred_to_be_copied = get(axes_to_be_copied,'children'); 
% Identify orientation of the axes 
[az,el] = view; 
% Copy the children of the axes 
copyobj(chilred_to_be_copied,a1); 
% Set the limits and orientation of the subplot as the original figure 
set(a1,'Xlim',get(axes_to_be_copied,'XLim')) 
set(a1,'Ylim',get(axes_to_be_copied,'YLim')) 
set(a1,'Zlim',get(axes_to_be_copied,'ZLim')) 
% set(a1,'title',title);
% title('u = 0.000');

% title(f_c.Title.String)
% One may set other properties such as labels, ticks etc. using the same 
% idea 
% Close the figure 
close(f_c); 
%% Copy second existing figure on second subplot 
% This figure contains a legend whose type is 'axes' hence it would require
% special treatment
% Open exiting figure 
f_c = openfig('ExampleGraphs/Part2.fig'); 
% Identify axes to be copied 
axes_to_be_copied = findobj(f_c,'type','axes', '-not', 'tag', 'legend');
% Identify the Legend
legend_to_be_copied = findobj(f_c,'type','axes', 'tag', 'legend');
% Identify the children of this axes 
chilred_to_be_copied = get(axes_to_be_copied,'children'); 
% Identify orientation of the axes 
[az,el] = view; 
% Copy the children of the axes 
copyobj(chilred_to_be_copied,a2);
% If there is a legend
if isfloat(legend_to_be_copied)
    copyobj(legend_to_be_copied, f);
end
% Set the limits and orientation of the subplot as the original figure 
set(a2,'Xlim',get(axes_to_be_copied,'XLim')) 
set(a2,'Ylim',get(axes_to_be_copied,'YLim')) 
set(a2,'Zlim',get(axes_to_be_copied,'ZLim')) 
view(a2,[az,el]) 
% One may set other properties such as labels, ticks etc. using the same 
% idea 
% Close the figure 
close(f_c); 
%% Copy second existing figure on second subplot 
% This figure contains a legend whose type is 'axes' hence it would require
% special treatment
% Open exiting figure 
f_c = openfig('ExampleGraphs/Part3.fig'); 
% Identify axes to be copied 
axes_to_be_copied = findobj(f_c,'type','axes', '-not', 'tag', 'legend');
% Identify the Legend
legend_to_be_copied = findobj(f_c,'type','axes', 'tag', 'legend');
% Identify the children of this axes 
chilred_to_be_copied = get(axes_to_be_copied,'children'); 
% Identify orientation of the axes 
[az,el] = view; 
% Copy the children of the axes 
copyobj(chilred_to_be_copied,a3);
% If there is a legend
if isfloat(legend_to_be_copied)
    copyobj(legend_to_be_copied, f);
end
% Set the limits and orientation of the subplot as the original figure 
set(a3,'Xlim',get(axes_to_be_copied,'XLim')) 
set(a3,'Ylim',get(axes_to_be_copied,'YLim')) 
set(a3,'Zlim',get(axes_to_be_copied,'ZLim')) 
view(a3,[az,el]) 
% One may set other properties such as labels, ticks etc. using the same 
% idea 
% Close the figure 
close(f_c); 
%% Copy second existing figure on second subplot 
% This figure contains a legend whose type is 'axes' hence it would require
% special treatment
% Open exiting figure 
f_c = openfig('ExampleGraphs/Part4.fig'); 
% Identify axes to be copied 
axes_to_be_copied = findobj(f_c,'type','axes', '-not', 'tag', 'legend');
% Identify the Legend
legend_to_be_copied = findobj(f_c,'type','axes', 'tag', 'legend');
% Identify the children of this axes 
chilred_to_be_copied = get(axes_to_be_copied,'children'); 
% Identify orientation of the axes 
[az,el] = view; 
% Copy the children of the axes 
copyobj(chilred_to_be_copied,a4);
% If there is a legend
if isfloat(legend_to_be_copied)
    copyobj(legend_to_be_copied, f);
end
% Set the limits and orientation of the subplot as the original figure 
set(a4,'Xlim',get(axes_to_be_copied,'XLim')) 
set(a4,'Ylim',get(axes_to_be_copied,'YLim')) 
set(a4,'Zlim',get(axes_to_be_copied,'ZLim')) 

view(a4,[az,el]) 
% One may set other properties such as labels, ticks etc. using the same 
% idea 
% Close the figure 
close(f_c); 