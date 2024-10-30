clear;
close all;
clc;

points = readtable('points.csv');
points1 = readtable('points1.csv');


list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'none');
end


figure(1)
plot(points.Var2(1), points.Var1(1),'o')
hold on
for i = 2:length(points.Var1)

    plot(points.Var2(i), points.Var1(i),'o')

end
legend('sa', 's1', 's2', 's3', 'e1', 'e2', 'e3', 'w1', 'w2', 'f1', 'f2', 'f3')

wing_outline_y = [points.Var1(1),  points.Var1(8), points.Var1(12), points.Var1(11), points.Var1(10), -0.08];
wing_outline_x = [points.Var2(1),  points.Var2(8), points.Var2(12), points.Var2(11), points.Var2(10), 0];



figure(2)
plot([points.Var2(3), points.Var2(7)], [points.Var1(3), points.Var1(7)]     ,'LineWidth',2)
hold on
plot([points.Var2(2), points.Var2(4)], [points.Var1(2), points.Var1(4)]     ,'LineWidth',2)
plot([points.Var2(1), points.Var2(6)], [points.Var1(1), points.Var1(6)]     ,'LineWidth',2)
plot([points.Var2(5), points.Var2(8)], [points.Var1(5), points.Var1(8)]     ,'LineWidth',2)
plot([points.Var2(7), points.Var2(9)], [points.Var1(7), points.Var1(9)]     ,'LineWidth',2)
plot([points.Var2(8), points.Var2(12)], [points.Var1(8), points.Var1(12)]   ,'LineWidth',2)
plot([points.Var2(8), points.Var2(10)], [points.Var1(8), points.Var1(10)]   ,'LineWidth',2)
plot([points.Var2(8), points.Var2(11)], [points.Var1(8), points.Var1(11)]   ,'LineWidth',2)

plot(wing_outline_x, wing_outline_y, 'LineStyle', '--','LineWidth',2, 'Color', '#808080')
xlabel('Span (m)')
ylabel('Chord (m)')
fontname(gca,"Times New Roman")
fontsize(gca, 12,"points")
xlim([0,0.6])
ylim([-0.2,0.1])
daspect([1 1 1])
legend('I','II','III','IV','V','VI','FI','FII')

figure(3)
plot(wing_outline_x, wing_outline_y,'LineStyle', '--', 'LineWidth',2, 'Color', '#808080')
hold on

xp1 = [points.Var2(8)/6:points.Var2(8)/6:points.Var2(8), points.Var2(8)+(points.Var2(12) - points.Var2(8))/6:(points.Var2(12) - points.Var2(8))/6:points.Var2(12)]            ;
yp1le = interp1([points.Var2(1),  points.Var2(8), points.Var2(12)],...
                [points.Var1(1),  points.Var1(8), points.Var1(12)],...
                xp1);

yp1te = interp1([0.0  , points.Var2(10), points.Var2(11), points.Var2(12)],...
                [-0.15, points.Var1(10), points.Var1(11), points.Var1(12)],...
                xp1);



for i = 1:12
    plot([xp1(i) xp1(i)], [yp1le(i) yp1te(i)], 'Color', '#000000', 'LineWidth',2)
end


xlabel('Span (m)')
ylabel('Chord (m)')
fontname(gca,"Times New Roman")
fontsize(gca, 12,"points")
xlim([0,0.6])
ylim([-0.3,0.1])
daspect([1 1 1])










    
wing_outline_y = [points1.Var1(1),  points1.Var1(8), points1.Var1(9), points1.Var1(12), points1.Var1(11), points1.Var1(10), -0.15];
wing_outline_x = [points1.Var2(1),  points1.Var2(8), points1.Var2(9), points1.Var2(12), points1.Var2(11), points1.Var2(10), 0];

figure(4)
plot([points1.Var2(1), points1.Var2(6)], [points1.Var1(1), points1.Var1(6)]     ,'LineWidth',2)
hold on
plot([points1.Var2(2), points1.Var2(4)], [points1.Var1(2), points1.Var1(4)]     ,'LineWidth',2)
plot([points1.Var2(3), points1.Var2(7)], [points1.Var1(3), points1.Var1(7)]     ,'LineWidth',2)
plot([points1.Var2(5), points1.Var2(8)], [points1.Var1(5), points1.Var1(8)]     ,'LineWidth',2)
plot([points1.Var2(7), points1.Var2(9)], [points1.Var1(7), points1.Var1(9)]     ,'LineWidth',2)
plot([points1.Var2(8), points1.Var2(10)], [points1.Var1(8), points1.Var1(10)]   ,'LineWidth',2)
plot([points1.Var2(8), points1.Var2(11)], [points1.Var1(8), points1.Var1(11)]   ,'LineWidth',2)
plot([points1.Var2(8), points1.Var2(12)], [points1.Var1(8), points1.Var1(12)]   ,'LineWidth',2)
plot(wing_outline_x, wing_outline_y,'LineStyle', '--', 'LineWidth',2, 'Color', '#808080')
xlabel('Span (m)')
ylabel('Chord (m)')
fontname(gca,"Times New Roman")
fontsize(gca, 12,"points")
daspect([1 1 1])
legend('1','2','3','4','5','6','7','8')


figure(5)
plot(wing_outline_x, wing_outline_y,'LineStyle', '--', 'LineWidth',2, 'Color', '#808080')
hold on

xp1 = [points1.Var2(8)/6:points1.Var2(8)/6:points1.Var2(8), points1.Var2(8)+(points1.Var2(12) - points1.Var2(8))/6:(points1.Var2(12) - points1.Var2(8))/6:points1.Var2(12)]            ;
yp1le = interp1([points1.Var2(1),  points1.Var2(8), points1.Var2(12)],...
                [points1.Var1(1),  points1.Var1(8), points1.Var1(12)],...
                xp1);

yp1te = interp1([0.0  , points1.Var2(10), points1.Var2(11), points1.Var2(12)],...
                [-0.15, points1.Var1(10), points1.Var1(11), points1.Var1(12)],...
                xp1);



for i = 1:12
    plot([xp1(i) xp1(i)], [yp1le(i) yp1te(i)], 'Color', '#000000', 'LineWidth',2)
end


xlabel('Span (m)')
ylabel('Chord (m)')
fontname(gca,"Times New Roman")
fontsize(gca, 12,"points")
daspect([1 1 1])


