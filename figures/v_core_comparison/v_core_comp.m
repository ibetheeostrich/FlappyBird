clear;
close all;
clc;

v_core = 0.0;

gamma = 0.2;

x = linspace(-1,1,400);
y = linspace(-1,1,400);

[X,Y] = meshgrid(x,y);

uind =  0.5 * pi*gamma * (Y) ./ sqrt(((X.^2) + (Y.^2)).^2 + v_core^4);
vind = -0.5 * pi*gamma * (X) ./ sqrt(((X.^2) + (Y.^2)).^2 + v_core^4);

Vind = sqrt(uind.^2 + vind.^2);

figure;
s = surf(X,Y,Vind);
s.EdgeColor = 'none';
c = colorbar('southoutside');
c.Label.String = 'Velocity';

axis equal
fontname(gcf,"Times New Roman")
fontsize(gcf,12,'points')
set(gca,'position',[0.2 0.3 0.7 0.75])
set(gcf,'units','centimeters','position',[10,10,7,8.5])
view(0,90)

xlabel('x')
ylabel('y')
