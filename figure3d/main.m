clearvars
tic
[Sample,Stochastic] = Parameter();
TeffRL = zeros( 1, Stochastic.aver );
Voltage_real = zeros( 6, Stochastic.aver );
energy_current = zeros( Sample.NWid * Sample.NLen, Stochastic.aver );
T_i = zeros( Sample.NWid * Sample.NLen, Stochastic.aver );
F_ii_ee = zeros( Sample.distribute_num_site, length( Sample.distribute_EF ), Stochastic.aver );
stream_seed = sum( clock );
for ii = 1:Stochastic.aver
    stream = RandStream( 'mlfg6331_64', 'Seed', stream_seed );
    RandStream.setGlobalStream( stream );
    stream.Substream = ii;
    [TeffRL( 1, ii ),Voltage_real( :, ii ),energy_current( :, ii ),T_i( :, ii ),F_ii_ee( :, :, ii )] = dissipation();
end
energy_current_aver = mean( energy_current, 2 );
T_i_aver = mean( T_i, 2 );
F_ii_ee_aver = mean( F_ii_ee, 3 );
save dissipation_data.mat TeffRL Voltage_real energy_current_aver T_i_aver F_ii_ee_aver
toc
