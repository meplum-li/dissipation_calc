%%% 边缘势、散射势垒、耗散源（只随机位置，固定耗散强度）
clearvars
tic
%以下参数的说明，可以直接参看Parameter()函数内的注释
[Sample, Stochastic] = Parameter();
parpool(Stochastic.core);
%%记录结果
aver = Stochastic.aver;
n_EF = length(Sample.hall_EF);
% TeffRL=zeros(length(Sample.hall_EF),Stochastic.aver);%两端的等效电导
TeffRL=repmat(zeros(n_EF,1),aver,1);%两端的等效电导
Voltage_real = repmat(zeros(6,n_EF),1,aver);%六个真实电极的电压

EF = repmat(Sample.hall_EF(:), aver, 1);
stream_seed = sum(clock);%随机数的设置
parfor ii = 1 : Stochastic.aver * length(Sample.hall_EF)
    stream=RandStream('mlfg6331_64', 'Seed', stream_seed);%随机数的设置
    RandStream.setGlobalStream(stream);%随机数的设置
    stream.Substream=ii;%随机数的设置
    [TeffRL(ii), Voltage_real(:, ii)] = dissipation(EF(ii));
end

Voltage_real = reshape(Voltage_real,[6, n_EF, aver]);%6 by energies by stochastic aver
TeffRL = reshape(TeffRL,[n_EF, aver]);%energies by stochastic aver
save dissipation_data.mat TeffRL Voltage_real
toc