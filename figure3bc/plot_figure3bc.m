clearvars
tic
load dissipation_data.mat
[Sample,Stochastic] = Parameter();
Lattice = hexagonal_lattice( Sample.Lx, Sample.Ly, 0 );
all_sites = sortrows( [ Lattice.siteA; Lattice.siteB ] );
energy_current_aver( all_sites( :, 5 )>0 & all_sites( :, 5 )<7 ) = 0;
energy_current_aver = reshape( energy_current_aver, Sample.NWid, Sample.NLen );
energy_current_x = reshape( all_sites( :, 1 ), [  ], Sample.NLen );
energy_current_y = reshape( all_sites( :, 2 ), [  ], Sample.NLen );
figure
surf( energy_current_x, energy_current_y, energy_current_aver )
shading interp ;
colormap jet ;
view( 0, 90 );
colorbar;
axis equal
xlim( Sample.Lx + [ -2, 2 ] )
ylim( Sample.Ly + [ -2, 2 ] )
clim( [ 0, 1E-6 ] )
set( gca, 'FontSize', 20 );
print( 'dissipated_energy_current', '-dpng', '-r200' )
distribut_site = sortrows( all_sites( all_sites( :, 6 )>0, : ) );
EE = Sample.distribute_EF;
cm = colormap( jet( 20 ) );
t = zeros( 10, 1 );
for ii = 2:2:20
    plot( EE, F_ii_ee_aver( ii, : ), "Color", cm( ii, : ), 'LineWidth', 3, 'DisplayName', [ '(x,y)=(', num2str( distribut_site( ii, 1 ) ), ', ', num2str( distribut_site( ii, 2 ) ), ')' ] )
    hold on
end
ylim( [ -0.1, 1.1 ] )
xlim( [ min( EE ), max( EE ) ] )
xlabel( 'E' )
ylabel( 'distribution' )
set( gca, 'FontSize', 20 );
legend( 'Show', 'FontSize', 16, 'Location', 'best' )
print( 'distribution', '-dpng', '-r200' )
