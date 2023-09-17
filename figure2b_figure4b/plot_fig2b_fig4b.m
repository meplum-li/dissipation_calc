clearvars
tic
load dissipation_data.mat
[Sample, Stochastic] = Parameter();
Lattice = hexagonal_lattice(Sample.Lx,Sample.Ly,0);
%%
%%%plot figure 2b
all_sites = sortrows([Lattice.siteA;Lattice.siteB]);
energy_current_aver = reshape(energy_current_aver,Sample.NWid,Sample.NLen);
energy_current_x = reshape(all_sites(:,1),[],Sample.NLen);
energy_current_y = reshape(all_sites(:,2),[],Sample.NLen);
figure
surf(energy_current_x,energy_current_y,energy_current_aver)
shading interp;
colormap jet;
view(0,90);
colorbar;
axis equal
xlim(Sample.Lx+[-2,2])
ylim(Sample.Ly+[-2,2])
clim([0, 1.9E-5])
set(gca, 'FontSize', 20);
print('dissipated_energy_current','-dpng','-r200')
%%
%%% plot figure 4b
all_sites = sortrows([Lattice.siteA;Lattice.siteB]);
T_i_aver = reshape(T_i_aver,Sample.NWid,Sample.NLen);
T_i_aver = T_i_aver*3.6*1000;%include the bias V_2p=3.6mV, and convert \Delta T to micro-Kelvin by multiplying 1000;
T_i_x = reshape(all_sites(:,1),[],Sample.NLen);
T_i_y = reshape(all_sites(:,2),[],Sample.NLen);
figure
surf(T_i_x,T_i_y,T_i_aver)
shading interp;
colormap jet;
view(0,90);
c=colorbar;
c.Label.String = '\DeltaT(mK)';
axis equal
xlim(Sample.Lx+[-2,2])
ylim(Sample.Ly+[-2,2])
clim([0, 4])
set(gca, 'FontSize', 20);
print('temperature_of_dissipation_source','-dpng','-r200')