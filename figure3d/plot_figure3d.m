clearvars
tic
load dissipation_data.mat
[Sample, Stochastic] = Parameter();
Lattice = hexagonal_lattice(Sample.Lx,Sample.Ly,0);
all_sites = sortrows([Lattice.siteA;Lattice.siteB]);

%% local distribution
distribut_site = sortrows(all_sites(all_sites(:,6)>0,:));
EE = Sample.distribute_EF;

cm = colormap(jet(20));        
t=zeros(10,1);% Define 0colormapt
% figure
for ii =1:1:10
    plot(EE, F_ii_ee_aver(ii,:), "Color", cm(ii,:), 'LineWidth',3,'DisplayName',['(x,y)=(',num2str(distribut_site(ii,1)),', ',num2str(distribut_site(ii,2)),')'])
    hold on
end
ylim([-0.1,1.1])
xlim([min(EE),max(EE)])
xlabel('E')
ylabel('distribution')
set(gca, 'FontSize', 20);
legend('Show','FontSize',16,'Location', 'best')
print('distribution','-dpng','-r200')