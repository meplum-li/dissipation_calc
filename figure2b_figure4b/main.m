clearvars
tic
[Sample, Stochastic] = Parameter();
parpool(Stochastic.core);
%%record results
TeffRL=zeros(1,Stochastic.aver);%effective conductance between the left and right contacts
Voltage_real = zeros(6,Stochastic.aver);%voltage of six real contacts
energy_current = zeros(Sample.NWid*Sample.NLen, Stochastic.aver);%spatial distribution of energy dissipation
T_i = zeros(Sample.NWid*Sample.NLen, Stochastic.aver);%spatial distribution of Buttiker probes
F_ii_ee = zeros(0, length(Sample.distribute_EF), Stochastic.aver);%energy distribution of electrons

stream_seed = sum(clock);%setting for random number
parfor ii = 1 : Stochastic.aver
    stream=RandStream('mlfg6331_64', 'Seed', stream_seed);%setting for random number
    RandStream.setGlobalStream(stream);%setting for random number
    stream.Substream=ii;%setting for random number
    [TeffRL(1, ii), Voltage_real(:, ii), energy_current(:, ii), T_i(:, ii), F_ii_ee(:,:,ii)] = dissipation();
end
energy_current_aver = mean(energy_current,2);
T_i_aver = mean(T_i,2);
F_ii_ee_aver = mean(F_ii_ee, 3);
save dissipation_data.mat TeffRL Voltage_real energy_current_aver T_i_aver F_ii_ee_aver
toc