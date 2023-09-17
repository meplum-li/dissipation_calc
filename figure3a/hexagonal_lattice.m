function Lattice = hexagonal_lattice(Lx,Ly,pltshow)
    [Sample,Stochastic] = Parameter();
    edge_wid = Sample.edge_wid;
    edge_wid_x = Sample.edge_wid_x;
    startA = [ 0, 0 ];
    dA = [ 0, 0 ];
    dB = [ 0, 1 / sqrt( 3 ) ];
    a1 = [ 1, 0 ];
    a2 = [ 1 / 2, sqrt( 3 ) / 2 ];
    bdyA = ([ a1', a2' ]) \ ([ kron( Lx, [ 1, 1 ] ); kron( [ 1, 1 ], Ly ) ] - kron( ones( 1, 4 ), startA' ));
    [ARi,ARj] = meshgrid( floor( min( bdyA( 1, : ) ) ):ceil( max( bdyA( 1, : ) ) ), floor( min( bdyA( 2, : ) ) ):ceil( max( bdyA( 2, : ) ) ) );
    siteA = startA + ARi( : ) * a1 + ARj( : ) * a2;
    idRA = (Lx( 1 ) - startA( 1 ) - siteA( :, 1 ))<=10 ^ -4 & (siteA( :, 1 ) - Lx( 2 ) + startA( 1 ))<=10 ^ -4 & (Ly( 1 ) - startA( 2 ) - siteA( :, 2 ))<=10 ^ -4 & (siteA( :, 2 ) - Ly( 2 ) + startA( 2 ))<=10 ^ -4;
    siteA = siteA( idRA, : );
    bdyB = ([ a1', a2' ]) \ ([ kron( Lx, [ 1, 1 ] ); kron( [ 1, 1 ], Ly ) ] - kron( ones( 1, 4 ), startA' + dB' ));
    [BRi,BRj] = meshgrid( floor( min( bdyB( 1, : ) ) ):ceil( max( bdyB( 1, : ) ) ), floor( min( bdyB( 2, : ) ) ):ceil( max( bdyB( 2, : ) ) ) );
    siteB = startA + dB + BRi( : ) * a1 + BRj( : ) * a2;
    idRB = (Lx( 1 ) - (startA( 1 )) - siteB( :, 1 ))<=10 ^ -3 & (siteB( :, 1 ) - Lx( 2 ) + (startA( 1 )))<=10 ^ -3 & (Ly( 1 ) - (startA( 2 )) - siteB( :, 2 ))<=10 ^ -3 & (siteB( :, 2 ) - Ly( 2 ) + (startA( 2 )))<=10 ^ -3;
    siteB = siteB( idRB, : );
    del1 = [ 0, 1 / sqrt( 3 ) ];
    del2 = [ 1 / 2, -1 / (2 * sqrt( 3 )) ];
    del3 = [ -1 / 2, -1 / (2 * sqrt( 3 )) ];
    tem_A = repelem( siteA, 3, 1 );
    A2B = tem_A + kron( ones( size( siteA, 1 ), 1 ), [ del1; del2; del3 ] );
    idA2B = ismembertol( A2B, siteB, 10 ^ -9, 'ByRows', true );
    bond = [ tem_A( idA2B, : ), A2B( idA2B, : ) ];
    Len = [ min( [ siteA( :, 1 ); siteB( :, 1 ) ] ), max( [ siteA( :, 1 ); siteB( :, 1 ) ] ) ];
    Wid = [ min( [ siteA( :, 2 ); siteB( :, 2 ) ] ), max( [ siteA( :, 2 ); siteB( :, 2 ) ] ) ];
    if pltshow
        figure( 'Visible', pltshow, 'Position', [ 0, 0, (Len( 2 ) - Len( 1 )) * 50, (Wid( 2 ) - Wid( 1 )) * 50 ] )
        cc = rand( 1, 3 );
        plot( siteA( :, 1 ), siteA( :, 2 ), 'o', 'MarkerSize', 4, 'MarkerFaceColor', cc, 'MarkerEdgeColor', cc, 'DisplayName', 'Sublattice A' );
        hold on
        cc = abs( 1 - cc );
        plot( siteB( :, 1 ), siteB( :, 2 ), 'o', 'MarkerSize', 4, 'MarkerFaceColor', cc, 'MarkerEdgeColor', cc, 'DisplayName', 'Sublattice B' );
        legend( 'AutoUpdate', 'off' )
        for ii = 1:size( bond, 1 )
            hold on
            plot( bond( ii, 1:2:3 ), bond( ii, 2:2:4 ), '-k', 'LineWidth', 1 )
        end
        xlim( Len )
        ylim( Wid )
        daspect( [ 1, 1, 1 ] )
    end
    Lattice.siteA = zeros( size( siteA, 1 ), 6 );
    Lattice.siteB = zeros( size( siteB, 1 ), 6 );
    Lattice.siteA( :, 1:3 ) = [ siteA, ones( size( siteA, 1 ), 1 ) ];
    Lattice.siteB( :, 1:3 ) = [ siteB, -ones( size( siteB, 1 ), 1 ) ];
    real_Wid = Wid( 2 ) - Wid( 1 );
    real_Len = Len( 2 ) - Len( 1 );
    Lattice.siteA( :, 4 ) = ~((Ly( 1 ) + edge_wid * real_Wid<=Lattice.siteA( :, 2 )) & (Lattice.siteA( :, 2 )<=Ly( 2 ) - edge_wid * real_Wid) & (Lx( 1 ) + edge_wid_x * real_Len<=Lattice.siteA( :, 1 )));
    Lattice.siteB( :, 4 ) = ~((Ly( 1 ) + edge_wid * real_Wid<=Lattice.siteB( :, 2 )) & (Lattice.siteB( :, 2 )<=Ly( 2 ) - edge_wid * real_Wid) & (Lx( 1 ) + edge_wid_x * real_Len<=Lattice.siteB( :, 1 )));
    topsl_x = Sample.topsl_x;
    topsl_y = Sample.topsl_y;
    botsl_x = Sample.botsl_x;
    botsl_y = Sample.botsl_y;
    sl_wid = Sample.sl_wid;
    Lattice.siteA( :, 5 ) = Lattice.siteA( :, 5 ) + (Lattice.siteA( :, 1 )==Lx( 1 ));
    Lattice.siteB( :, 5 ) = Lattice.siteB( :, 5 ) + (Lattice.siteB( :, 1 )==Lx( 1 ));
    Lattice.siteA( :, 5 ) = Lattice.siteA( :, 5 ) + 4 * (Lattice.siteA( :, 1 )==Lx( 2 ));
    Lattice.siteB( :, 5 ) = Lattice.siteB( :, 5 ) + 4 * (Lattice.siteB( :, 1 )==Lx( 2 ));
    A_possible_dissipation = find( Lattice.siteA( :, 5 )==0 );
    id_A_dissipation = randperm( length( A_possible_dissipation ), round( length( A_possible_dissipation ) * Stochastic.dissipation.dens ) );
    Lattice.siteA( A_possible_dissipation( id_A_dissipation ), 5 ) = 7;
    B_possible_dissipation = find( Lattice.siteB( :, 5 )==0 );
    id_B_dissipation = randperm( length( B_possible_dissipation ), round( length( B_possible_dissipation ) * Stochastic.dissipation.dens ) );
    Lattice.siteB( B_possible_dissipation( id_B_dissipation ), 5 ) = 7;
    Lattice.num_dissipation = nnz( Lattice.siteA( :, 5 )==7 ) + nnz( Lattice.siteB( :, 5 )==7 );
    Lattice.bond = bond;
    Lattice.Len = Len;
    Lattice.Wid = Wid;
    fprintf( "length of the graphene is %.2fa\n", Len( 2 ) - Len( 1 ) )
    fprintf( "width of the graphene is %.2fa\n", Wid( 2 ) - Wid( 1 ) )
end
