clearvars
tic
[Sample,Stochastic] = Parameter();
parpool( Stochastic.core );
aver = Stochastic.aver;
n_EF = length( Sample.hall_EF );
TeffRL = repmat( zeros( n_EF, 1 ), aver, 1 );
Voltage_real = repmat( zeros( 6, n_EF ), 1, aver );
EF = repmat( Sample.hall_EF( : ), aver, 1 );
stream_seed = sum( clock );
parfor ii = 1:Stochastic.aver * length( Sample.hall_EF )
    stream = RandStream( 'mlfg6331_64', 'Seed', stream_seed );
    RandStream.setGlobalStream( stream );
    stream.Substream = ii;
    [TeffRL( ii ),Voltage_real( :, ii )] = dissipation( EF( ii ) );
end
Voltage_real = reshape( Voltage_real, [ 6, n_EF, aver ] );
TeffRL = reshape( TeffRL, [ n_EF, aver ] );
save dissipation_data.mat TeffRL Voltage_real
toc
