figure(1)
plot(test1.e2,test1.e00,'LineWidth',2)
hold on
plot(wttest1.Times,wttest1.DragN,'LineWidth',2)
hold off
set(gca,'fontname','TW Cen MT') 
set(gca,'FontSize',16)
grid on
grid minor
xlim([0,0.33])
xlabel('Time (s)')
ylabel('Drag (N)')
legend('LDVM','Wind Tunnel')

figure(2)
plot(test2.e1(4:end),test2.e01(4:end),'LineWidth',2)
hold on
plot(wttest2.Times,wttest2.DragN,'LineWidth',2)
hold off
set(gca,'fontname','TW Cen MT') 
set(gca,'FontSize',16)
grid on
grid minor
xlabel('Time (s)')
ylabel('Drag (N)')
legend('LDVM','Wind Tunnel')