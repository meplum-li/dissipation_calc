clearvars
tic
load dissipation_data.mat
[Sample, Stochastic] = Parameter();
Lattice = hexagonal_lattice(Sample.Lx,Sample.Ly,0);
%%
R23 = (Voltage_real(2,:) - Voltage_real(3,:))./(TeffRL.');
R26 =  (Voltage_real(2,:) - Voltage_real(6,:))./(TeffRL.');
figure
yyaxis right
plot(Sample.hall_EF(:), R23(:),'Color',rand(1,3),'LineWidth',2,'DisplayName','R_{23}')
ylabel('R_{xx}')
hold on
yyaxis left
plot(Sample.hall_EF(:), R26(:),'Color',rand(1,3),'LineWidth',2,'DisplayName','R_{26}')
hold on
xlim([min(Sample.hall_EF(:)),Sample.hall_EF(end)])
ylim([-1.1,1.1])
xlabel('E_F')
ylabel('R_{H}')
set(gca, 'FontSize', 15);
legend('Show','FontSize',12,'Location', 'best')
print('Hall_measurement','-dpng','-r200')
%%
R2p = (Voltage_real(1,:) - Voltage_real(4,:))./(TeffRL.');
figure
% yyaxis left
plot(Sample.hall_EF(:), R2p(:),'Color',rand(1,3),'LineWidth',2,'DisplayName','R_{2p}')
ylabel('R_{2p}')
hold on
xlim([min(Sample.hall_EF(:)),Sample.hall_EF(end)])
ylim([0,1.1])
xlabel('E_F')
ylabel('R_{2p}')
set(gca, 'FontSize', 15);
legend('Show','FontSize',16,'Location', 'best')
print('two-terminal_resistance','-dpng','-r200')