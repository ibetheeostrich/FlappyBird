clear;
close all;
clc;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'none');
end
set(groot,...
'defaulttextinterpreter','latex',...
'defaultAxesTickLabelInterpreter','latex',...
'defaultLegendInterpreter','latex')

file = '0.3000_5.0deg_lift.csv';
name = file(1:end-4);


imp = readmatrix(file);

u = [8 9 10 11 12 13 14];
hz = [1 1.5 2 2.5 3];

[U, H] = meshgrid(u,hz);

figure;
s = contourf(U,H,imp.*H,'LineStyle','none');
c = colorbar('southoutside');
c.Label.String = 'Lift Impulse ($$Ns$$)';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

% fontname(gcf,"Times New Roman")
fontsize(gcf,12,'points')
% set(gca,'position',[0.2 0.3 0.7 0.75])
xlabel("Velocity ($$m/s$$)",'FontSize',12)
ylabel("Wing Beat Frequency (Hz)",'FontSize',12)
set(gcf,'units','centimeters','position',[10,10,7.5,10])
xlim([8,14])
ylim([1,3])
% view(0,90)
ax = gca;
% exportgraphics(ax,strcat(name,'isecimp','.pdf'),'ContentType','vector')

