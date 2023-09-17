function nearest = isnearest(uc,nuc)
    [xa,xb] = meshgrid( nuc( :, 1 ), uc( :, 1 ) );
    delx = xa - xb;
    [ya,yb] = meshgrid( nuc( :, 2 ), uc( :, 2 ) );
    dely = ya - yb;
    tol = 10 ^ -4;
    nearest.id = abs( sqrt( delx .^ 2 + dely .^ 2 ) - 1 / sqrt( 3 ) )<tol;
    nearest.xa = xa;
    nearest.xb = xb;
    nearest.ya = ya;
    nearest.yb = yb;
end
