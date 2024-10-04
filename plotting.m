clear;
close all;
clc;

points = readtable('points.csv');



figure(1)
plot(points.Var2(1), points.Var1(1),'o')
hold on
for i = 2:length(points.Var1)

    plot(points.Var2(i), points.Var1(i),'o')

end
legend('sa', 's1', 's2', 's3', 'e1', 'e2', 'e3', 'w1', 'w2', 'f1', 'f2', 'f3')

figure(2)
plot([points.Var2(1), points.Var2(6)], [points.Var1(1), points.Var1(6)])
hold on
plot([points.Var2(2), points.Var2(4)], [points.Var1(2), points.Var1(4)])
plot([points.Var2(3), points.Var2(7)], [points.Var1(3), points.Var1(7)])
plot([points.Var2(5), points.Var2(8)], [points.Var1(5), points.Var1(8)])
plot([points.Var2(7), points.Var2(9)], [points.Var1(7), points.Var1(9)])
plot([points.Var2(8), points.Var2(10)], [points.Var1(8), points.Var1(10)])
plot([points.Var2(8), points.Var2(11)], [points.Var1(8), points.Var1(11)])
plot([points.Var2(8), points.Var2(12)], [points.Var1(8), points.Var1(12)])

