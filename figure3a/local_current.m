tic
[Sample,Stochastic] = Parameter();
Lattice = hexagonal_lattice( Sample.Lx, Sample.Ly, 0 );
pltshow = 0;
if pltshow
    figure
    all_sites = [ Lattice.siteA; Lattice.siteB ];
    for ii = 1:6
        plot_sites = all_sites( all_sites( :, 5 )==ii, 1:2 );
        scatter( plot_sites( :, 1 ), plot_sites( :, 2 ), 10, rand( 1, 3 ), 'o', 'filled', 'MarkerFaceAlpha', 1 )
        hold on
    end
    plot_sites = all_sites( all_sites( :, 5 )==7, 1:2 );
    scatter( plot_sites( :, 1 ), plot_sites( :, 2 ), 20, 'k', 'x' )
    hold on
    plot_sites = all_sites( all_sites( :, 6 )==1, 1:2 );
    scatter( plot_sites( :, 1 ), plot_sites( :, 2 ), 20, 'k', 's' )
    xlim( Sample.Lx + [ -2, 2 ] )
    ylim( Sample.Ly + [ -2, 2 ] )
    set( gca, 'FontSize', 20 );
    daspect( [ 1, 1, 1 ] )
    print( 'device', '-dpng', '-r300' )
end
[GAMMAlead_left,GAMMAlead_right,SIGMAlead_left,SIGMAlead_right] = lead( Sample.EF, Sample.phiz, Lattice );
SIGMAsl = -1i / 2 * Sample.GAMMAsl;
siteA = Lattice.siteA;
siteB = Lattice.siteB;
sa = min( siteA( :, 1 ) );
edgeA = siteA( siteA( :, 1 )==sa, : );
edgeB = siteB( siteB( :, 1 )==sa, : );
uc = sortrows( [ edgeA; edgeB ], 2 );
res_sites = [ Lattice.siteA; Lattice.siteB ];
N = size( res_sites, 1 );
nu = {  };
nu{ 1 } = uc;
temp_N = size( nu{ 1 }, 1 );
N_nu = [  ];
N_nu( 1 ) = temp_N;
ii = 1;
res_sites = res_sites( ~ismembertol( res_sites, nu{ 1 }, 10 ^ -9, 'ByRows', true ), : );
while temp_N<N
    last_nu = nu{ ii };
    [temp_nu,res_sites] = nextslice( last_nu, res_sites );
    nu{ ii + 1 } = temp_nu;
    ii = ii + 1;
    N_nu( ii ) = size( temp_nu, 1 );
    temp_N = temp_N + N_nu( ii );
end
Nslice = size( nu, 2 );
pltshow = 0;
if pltshow
    figure( 'Visible', 1, 'Position', [ 0, 0, (Sample.Lx( 2 ) - Sample.Lx( 1 )) * 50, (Sample.Ly( 2 ) - Sample.Ly( 1 )) * 50 ] )
    for ii = 1:Nslice
        scatter( nu{ ii }( :, 1 ), nu{ ii }( :, 2 ), 10, rand( 1, 3 ), 'o', 'filled', 'MarkerFaceAlpha', 1 )
        hold on
        scatter( nu{ ii }( logical( nu{ ii }( :, 4 ) ), 1 ), nu{ ii }( logical( nu{ ii }( :, 4 ) ), 2 ), 16, 'r', 's', 'filled', 'MarkerFaceAlpha', 0.5 )
        hold on
    end
    for ii = 1:size( Lattice.bond, 1 )
        plot( Lattice.bond( ii, 1:2:3 ), Lattice.bond( ii, 2:2:4 ), '-k', 'LineWidth', 1 )
        hold on
    end
    xlim( Sample.Lx + [ -2, 2 ] )
    ylim( Sample.Ly + [ -2, 2 ] )
    set( gca, 'FontSize', 20 );
    daspect( [ 1, 1, 1 ] )
    print( 'device', '-dpng', '-r100' )
end
H00 = zeros( N_nu( 1 ), N_nu( 1 ), Nslice );
H01 = cell( 2, 1 );
nearest00 = cell( 2, 1 );
nearest00{ 1 } = isnearest( nu{ 1 }, nu{ 1 } );
nearest01_1 = isnearest( nu{ 1 }, nu{ 2 } );
H01{ 1 } = sparse( nearest01_1.id * (-Sample.t) .* exp( 1i * 2 * pi * Sample.phiz * (nearest01_1.xa - nearest01_1.xb) .* (nearest01_1.ya + nearest01_1.yb) / 2 ) );
nearest00{ 2 } = isnearest( nu{ 2 }, nu{ 2 } );
nearest01_2 = isnearest( nu{ 2 }, nu{ 3 } );
H01{ 2 } = sparse( nearest01_2.id * (-Sample.t) .* exp( 1i * 2 * pi * Sample.phiz * (nearest01_2.xa - nearest01_2.xb) .* (nearest01_2.ya + nearest01_2.yb) / 2 ) );
SIGMA_dissipation = zeros( Lattice.num_dissipation, 1 );
tem_num_diss = 0;
for ii = 1:1:Nslice
    nii = N_nu( ii );
    tem_tem_diss = nnz( nu{ ii }( :, 5 )==7 );
    tem_SIGMA_dissipation = Stochastic.dissipation.str * (nu{ ii }( :, 5 )==7);
    SIGMA_dissipation( tem_num_diss + 1:tem_num_diss + tem_tem_diss ) = tem_SIGMA_dissipation( nu{ ii }( :, 5 )==7 );
    tem_num_diss = tem_num_diss + tem_tem_diss;
    H00( :, :, ii ) = sparse( Sample.V0 * eye( nii ) + Sample.gamma1 * diag( nu{ ii }( :, 4 ) ) + Stochastic.disorder.str * spdiags( rand( nii, 1 ) - 0.5, 0, nii, nii ) + Sample.SQUID.str * spdiags( exp( -vecnorm( (nu{ ii }( :, 1:2 ) - Sample.SQUID.location), 2, 2 ) / Sample.SQUID.radius ), 0, nii, nii ) + (ii==1) * SIGMAlead_left + (ii==Nslice) * SIGMAlead_right + (ii~=1) * (ii~=Nslice) * SIGMAsl * diag( nu{ ii }( :, 5 )==2 | nu{ ii }( :, 5 )==3 | nu{ ii }( :, 5 )==5 | nu{ ii }( :, 5 )==6 ) + spdiags( tem_SIGMA_dissipation, 0, nii, nii ) + nearest00{ 2 - mod( ii, 2 ) }.id * (-Sample.t) .* exp( 1i * 2 * pi * Sample.phiz * (nearest00{ 2 - mod( ii, 2 ) }.xa - nearest00{ 2 - mod( ii, 2 ) }.xb) .* (nearest00{ 2 - mod( ii, 2 ) }.ya + nearest00{ 2 - mod( ii, 2 ) }.yb) / 2 ) );
end
SIGin = cell( Nslice, 1 );
SIGin{ Nslice } = zeros( N_nu( Nslice ) );
SIGi1 = cell( Nslice, 1 );
SIGi1{ 1 } = zeros( N_nu( 1 ) );
SIGi1_h = SIGi1;
for ii = Nslice - 1:-1:1
    SIGin{ ii } = H01{ -mod( ii, 2 ) + 2 } / (Sample.EF * speye( N_nu( ii + 1 ) ) - H00( :, :, ii + 1 ) - SIGin{ ii + 1 }) * H01{ -mod( ii, 2 ) + 2 }';
end
for ii = 2:Nslice
    SIGi1_h{ ii } = H01{ mod( ii, 2 ) + 1 }' / (Sample.EF * eye( N_nu( ii - 1 ) ) - H00( :, :, ii - 1 ) - SIGi1{ ii - 1 });
    SIGi1{ ii } = SIGi1_h{ ii } * H01{ mod( ii, 2 ) + 1 };
end
current_intra = zeros( N_nu( 1 ), N_nu( 1 ), Nslice );
current_inter = zeros( N_nu( 1 ), N_nu( 1 ), Nslice );
last_GRi1 = inv( Sample.EF * eye( N_nu( 1 ) ) - H00( :, :, 1 ) - SIGi1{ 1 } - SIGin{ 1 } );
prod_SIG_h = 1;
for ii = 1:Nslice - 1
    Gii_next = inv( Sample.EF * speye( N_nu( ii + 1 ) ) - H00( :, :, ii + 1 ) - SIGi1{ ii + 1 } - SIGin{ ii + 1 } );
    prod_SIG_h = SIGi1_h{ ii + 1 } * prod_SIG_h;
    new_GRi1 = Gii_next * prod_SIG_h;
    Hij = sparse( triu( nearest00{ 2 - mod( ii, 2 ) }.id ) ) * (-Sample.t) .* exp( 1i * 2 * pi * Sample.phiz * (nearest00{ 2 - mod( ii, 2 ) }.xa - nearest00{ 2 - mod( ii, 2 ) }.xb) .* (nearest00{ 2 - mod( ii, 2 ) }.ya + nearest00{ 2 - mod( ii, 2 ) }.yb) / 2 );
    current_intra( :, :, ii ) = 2 * imag( Hij .* (last_GRi1 * GAMMAlead_left * last_GRi1').' );
    current_inter( :, :, ii ) = 2 * imag( H01{ -mod( ii, 2 ) + 2 } .* (new_GRi1 * GAMMAlead_left * last_GRi1').' );
    last_GRi1 = new_GRi1;
end
ii = Nslice;
Hij = sparse( triu( nearest00{ 2 - mod( ii, 2 ) }.id ) ) * (-Sample.t) .* exp( 1i * 2 * pi * Sample.phiz * (nearest00{ 2 - mod( ii, 2 ) }.xa - nearest00{ 2 - mod( ii, 2 ) }.xb) .* (nearest00{ 2 - mod( ii, 2 ) }.ya + nearest00{ 2 - mod( ii, 2 ) }.yb) / 2 );
current_intra( :, :, ii ) = 2 * imag( Hij .* (last_GRi1 * GAMMAlead_left * last_GRi1').' );
save local_current.mat current_intra current_inter
toc
