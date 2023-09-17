function [Sample,Stochastic] = Parameter
    Sample.EF = 0.05;
    Sample.hall_EF = linspace( -0.9, 1.1, 32 * 7 );
    Sample.t = 2.7;
    Sample.V0 = 0;
    Sample.gamma2 = 0;
    Sample.M = 0;
    Sample.EFlead = 1.5 * Sample.t;
    Sample.eta = 10 ^ (-6);
    Sample.phiz = -0.004;
    Sample.B_real = -Sample.phiz * 2 * 2.0678 * 10 ^ (-15) / (0.246 * 10 ^ (-9)) ^ 2;
    Sample.Lx = [ 0, 202 + 1 / 2 ];
    Sample.Ly = [ -sqrt( 3 ) / 2 * 1 / 3, sqrt( 3 ) * (120 - 0.5) ];
    Sample.gamma1 = 0.1;
    Sample.edge_wid = 60 / (Sample.Ly( 2 ) - Sample.Ly( 1 ));
    Sample.edge_wid_x = 20 / (Sample.Lx( 2 ) - Sample.Lx( 1 ));
    Sample.NWid = 2 * round( (Sample.Ly( 2 ) + sqrt( 3 ) / 2) / sqrt( 3 ) );
    Sample.NLen = ((Sample.Lx( 2 ) - Sample.Lx( 1 ))) / (0.5) + 1;
    Sample.N = Sample.NLen * Sample.NWid;
    Sample.topsl_x = [ 30, 170 ] + 0.5;
    Sample.topsl_y = Sample.Ly( 2 );
    Sample.botsl_x = Sample.topsl_x;
    Sample.botsl_y = Sample.Ly( 1 );
    Sample.sl_wid = 9;
    Sample.GAMMAsl = 2 * pi * 2;
    Sample.SQUID.str = -0.13 * Sample.t;
    Sample.SQUID.radius = 10;
    Sample.SQUID.location = [ mean( Sample.topsl_x + [ Sample.sl_wid, 0 ] ), Sample.Ly( 1 ) + (Sample.Ly( 2 ) - Sample.Ly( 1 )) * (Sample.edge_wid / 2) ];
    Sample.distribute_EF = linspace( -0.1, 1.1, 100 );
    Stochastic.dissipation.str = -1i / 2 * 0.1 * Sample.t;
    Stochastic.dissipation.dens = 0.005;
    Stochastic.dissipation.temper = 0.1;
    Stochastic.disorder.str = 0 * Sample.t;
    Stochastic.disorder.dens = 1;
    Stochastic.aver = 50;
    Stochastic.core = 32;
end
