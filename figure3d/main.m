%%% 边缘势、散射势垒、耗散源（只随机位置，固定耗散强度）
clearvars
tic
%以下参数的说明，可以直接参看Parameter()函数内的注释
[Sample, Stochastic] = Parameter();
% parpool(Stochastic.core);
%%记录结果
TeffRL=zeros(1,Stochastic.aver);%两端的等效电导
Voltage_real = zeros(6,Stochastic.aver);%六个真实电极的电压
energy_current = zeros(Sample.NWid*Sample.NLen, Stochastic.aver);%能流的空间分布
T_i = zeros(Sample.NWid*Sample.NLen, Stochastic.aver);%虚拟电极温度的空间分布
F_ii_ee = zeros(Sample.distribute_num_site, length(Sample.distribute_EF), Stochastic.aver);%电子的分布

stream_seed = sum(clock);%随机数的设置
for ii = 1 : Stochastic.aver
    stream=RandStream('mlfg6331_64', 'Seed', stream_seed);%随机数的设置
    RandStream.setGlobalStream(stream);%随机数的设置
    stream.Substream=ii;%随机数的设置
    [TeffRL(1, ii), Voltage_real(:, ii), energy_current(:, ii), T_i(:, ii), F_ii_ee(:,:,ii)] = dissipation();
end
energy_current_aver = mean(energy_current,2);
T_i_aver = mean(T_i,2);
F_ii_ee_aver = mean(F_ii_ee, 3);
save dissipation_data.mat TeffRL Voltage_real energy_current_aver T_i_aver F_ii_ee_aver
toc