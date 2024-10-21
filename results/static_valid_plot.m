clear;
close all;
clc;

deg0_0 =    readtable('0deg.csv');
deg2_5 =    readtable('2.5deg.csv');
deg5_0 =    readtable('5deg.csv');
deg7_5 =    readtable('7.5deg.csv');
deg10 =     readtable('10deg.csv');


cl = [deg0_0.Var1,deg2_5.Var1,deg5_0.Var1,deg7_5.Var1,deg10.Var1];
td = deg0_0.Var2;

cl_static = [0, 2*pi*deg2rad(2.5), 2*pi*deg2rad(5.0), 2*pi*deg2rad(7.5), 2*pi*deg2rad(10)];

for i = 1:5

    figure(i)
    plot(td(4:end),cl(4:end,i), LineWidth=2,Color="b")
    hold on
    plot(td,zeros(length(td),1) + cl_static(i), "LineStyle","--",LineWidth=2,Color='k')

    p = patch([0 0  td(end) td(end) 0], [cl_static(i) 0.9*cl_static(i) 0.9*cl_static(i) cl_static(i) cl_static(i)],'r','FaceAlpha',0.2,'EdgeColor','none');
    ylim([0,1.5])
    grid on
    grid minor
    
    legend('LUATAT+DVM','Ideal Flat Plate')

    ylabel('Coefficient of Lift')
    xlabel('Time (s)')
end