%%% 计算六段体系的霍尔以及耗散效应
%随机耗散位置以及随机的耗散强度
tic
%=========%
%以下参数的说明，可以直接参看Parameter()函数内的注释
[Sample, Stochastic] = Parameter();
%%% 侧面电极的几何设置
%                |2|         |3|
%     ---===========---
%     -1-===========-4-
%     ---===========---
%                |6|         |5|
%%
%生成格子
Lattice = hexagonal_lattice(Sample.Lx,Sample.Ly,0);
%%% plot the partition scheme
pltshow = 0;
if pltshow
    figure
    all_sites=[Lattice.siteA;Lattice.siteB];
    for  ii = 1:6
        plot_sites = all_sites(all_sites(:,5)==ii, 1:2);
        scatter(plot_sites(:,1), plot_sites(:,2), 10, rand(1,3), 'o', 'filled', 'MarkerFaceAlpha',1)
        hold on
    end
    plot_sites = all_sites(all_sites(:,5)==7, 1:2);
    scatter(plot_sites(:,1), plot_sites(:,2), 20, 'k', 'x')
    hold on
    plot_sites = all_sites(all_sites(:,6)==1, 1:2);
    scatter(plot_sites(:,1), plot_sites(:,2), 20, 'k', 's')
    xlim(Sample.Lx+[-2,2])
    ylim(Sample.Ly+[-2,2])
    % axis tight;
    set(gca, 'FontSize', 20);
    daspect([1 1 1])
    print('device','-dpng','-r300')
end
%% 表面格林函数
%=========%左右导线的表面格林函数
[GAMMAlead_left,GAMMAlead_right, SIGMAlead_left, SIGMAlead_right]=lead(Sample.EF,Sample.phiz,Lattice);
%=========%
%=========%侧面电极的输运参数设置。侧面的电极参数均一致。
SIGMAsl = -1i/2 * Sample.GAMMAsl;
%=========%
%=========%虚拟导线的表面格林函数
%Parameter()函数内已经设置
%=========%

%% 生成first slice (只有side leads，通电流的lead是上下边中间的lead)
% [SA,LA]=bounds(Lattice.siteA,1);
% [SB,LB]=bounds(Lattice.siteB,1);
% Wid = Sample.Ly(2) - Sample.Ly(1);
siteA = Lattice.siteA;
siteB = Lattice.siteB;
sa = min(siteA(:,1));
% sb = min(Lattice.siteB(:,1));
edgeA = siteA(siteA(:,1)==sa,:);%最左边的A格点
edgeB = siteB(siteB(:,1)==sa,:);%最左边的B格点
uc = sortrows([edgeA;edgeB],2);%按照格点的y坐标，从小到大对最左边的元胞内的格点进行排序。uc的第三个分量标志子格子的类型，1 for A and -1 for B %unit cell
% uc也正好是左电极的位置
%% 划分迭代方案
res_sites=[Lattice.siteA;Lattice.siteB];%剩余未被分配的格点
N = size(res_sites,1);%总的格点数目
nu={};%初始化，格点的分配方案
nu{1}=uc;%第一层格点
temp_N = size(nu{1},1);%已经被分配的格点的数目
N_nu=[];% number of sites in each slice
N_nu(1)=temp_N;%当前层的格点数目
ii=1;%分配第一层
res_sites = res_sites(~ismembertol(res_sites,nu{1},10^-9,'ByRows',true),:);%剔除被分配到第一层的格点
while temp_N<N
    last_nu=nu{ii}; %记录上一层的格点
    [temp_nu, res_sites] = nextslice(last_nu, res_sites);
    nu{ii+1}=temp_nu;%记录第ii+1个slice的格点
    ii = ii+1;
    N_nu(ii) = size(temp_nu,1);%第ii+1层格点
    temp_N=temp_N+N_nu(ii);%更新已经被分配的格点的数目
end
Nslice = size(nu,2);
%% plot the partition scheme
pltshow = 0;
if pltshow
    figure('Visible',1,'Position',[0 0 (Sample.Lx(2)-Sample.Lx(1))*50 (Sample.Ly(2)-Sample.Ly(1))*50])
    for ii = 1 : Nslice
        scatter(nu{ii}(:,1),nu{ii}(:,2),10,rand(1,3),'o','filled','MarkerFaceAlpha',1)
        hold on
        scatter(nu{ii}(logical(nu{ii}(:,4)),1),nu{ii}(logical(nu{ii}(:,4)),2),16,'r','s','filled','MarkerFaceAlpha',0.5)
        hold on
    end
    for ii = 1 : size(Lattice.bond,1)
        plot(Lattice.bond(ii,1:2:3),Lattice.bond(ii,2:2:4),'-k','LineWidth',1)
        hold on
    end
    xlim(Sample.Lx+[-2,2])
    ylim(Sample.Ly+[-2,2])
    % axis tight;
    set(gca, 'FontSize', 20);
    daspect([1 1 1])
    print('device','-dpng','-r100')
end
%% 生成哈密顿量
%=======%生成哈密顿量
H00=zeros(N_nu(1),N_nu(1),Nslice);
H01 = cell(2,1);
nearest00 = cell(2, 1);
%这种特殊的情况下,slice内部只有两种不同的最近邻关系。由于边缘势垒影响在位能，因此需要记录每个slice的H00。
%slice之间的最近邻关系，也只有两种
%%% the first slice
nearest00{1} = isnearest(nu{1},nu{1});
nearest01_1 = isnearest(nu{1}, nu{2});%(nu{1})by(nu{2})
H01{1} = sparse(nearest01_1.id*(-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest01_1.xa - nearest01_1.xb).*(nearest01_1.ya + nearest01_1.yb)/2));
%%%the second kind of slice
nearest00{2}  = isnearest(nu{2},nu{2});
nearest01_2 = isnearest(nu{2}, nu{3});%(nu{2})by(nu{3})
H01{2} = sparse(nearest01_2.id*(-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest01_2.xa - nearest01_2.xb).*(nearest01_2.ya + nearest01_2.yb)/2));
SIGMA_dissipation = zeros(Lattice.num_dissipation,1);%记录随机的耗散源的自能
tem_num_diss = 0;
for ii = 1 : 1 : Nslice
    nii = N_nu(ii);
    tem_tem_diss = nnz(nu{ii}(:,5)==7);
    tem_SIGMA_dissipation = Stochastic.dissipation.str * ( nu{ii}(:,5)==7 );
    SIGMA_dissipation(tem_num_diss+1 : tem_num_diss + tem_tem_diss) = tem_SIGMA_dissipation(nu{ii}(:,5)==7);
    tem_num_diss = tem_num_diss + tem_tem_diss;
    H00(:,:,ii)= sparse( ...
        Sample.V0*eye(nii)...%在位能
        + Sample.gamma1*diag(nu{ii}(:,4))...%边缘势
        + Stochastic.disorder.str * spdiags(rand(nii,1) - 0.5, 0, nii, nii)...%安德森无序
        + Sample.SQUID.str * spdiags(exp(-vecnorm((nu{ii}(:,1:2) - Sample.SQUID.location),2,2)/Sample.SQUID.radius), 0, nii, nii)...%SQUID引起的potential %+ Sample.SQUID.str * spdiags( prod(abs( nu{ii}(:,1:2)-Sample.SQUID.location )<[Sample.SQUID.Len/2, Sample.SQUID.Wid/2],2 ), 0, nii, nii)...%SQUID引起的potential
        + (ii == 1) * SIGMAlead_left + (ii == Nslice) * SIGMAlead_right...%左右电极
        + (ii ~= 1) * (ii ~= Nslice) * SIGMAsl * diag( nu{ii}(:,5)==2|nu{ii}(:,5)==3|nu{ii}(:,5)==5|nu{ii}(:,5)==6 )...%侧面电极
        + spdiags(tem_SIGMA_dissipation, 0, nii, nii)...%耗散电极
        + nearest00{2-mod(ii,2)}.id*(-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest00{2-mod(ii,2)}.xa - nearest00{2-mod(ii,2)}.xb).*(nearest00{2-mod(ii,2)}.ya + nearest00{2-mod(ii,2)}.yb)/2)...
        );
end
%% green function 这里只需要GR{i,i}，以及GR{i,1}
% SIGin=zeros(Nx,Nx,Ny); 只需要\Sigma_{1:Nslice}^{(Nslice)}
SIGin = cell(Nslice,1);
SIGin{Nslice}=zeros(N_nu(Nslice));
% SIGi1=zeros(Nx,Nx,Ny); 只需要\Sigma_{1:Nslice}^{(1)}
SIGi1=cell(Nslice,1);
SIGi1{1}=zeros(N_nu(1));
SIGi1_h = SIGi1;

for ii = Nslice-1:-1:1
    SIGin{ii}=H01{-mod(ii,2)+2}/( Sample.EF*speye(N_nu(ii+1))-H00(:,:,ii+1)-SIGin{ii+1})*H01{-mod(ii,2)+2}';
end
%\Sigma_i^{(1)}
% SIGi1_h{2} = H01{1}'*inv( EF*eye(length(nu{1})) -H00{1});
% SIGi1{2} = SIGi1_h{2}*H01{1};
for ii = 2 : Nslice
    SIGi1_h{ii} = H01{mod(ii,2)+1}'/( Sample.EF*eye(N_nu(ii-1)) - H00(:,:,ii-1)-SIGi1{ii-1});
    SIGi1{ii} = SIGi1_h{ii}*H01{mod(ii,2)+1};
end
% get GR and current
current_intra = zeros(N_nu(1),N_nu(1),Nslice);%层内电流
current_inter = zeros(N_nu(1),N_nu(1),Nslice);%层间电流
last_GRi1 = inv( Sample.EF*eye(N_nu(1)) -H00(:,:,1) - SIGi1{1} -SIGin{1});
prod_SIG_h = 1;%initialization
for ii = 1 : Nslice-1
    Gii_next = inv( Sample.EF*speye(N_nu(ii+1)) - H00(:,:,ii+1) - SIGi1{ii+1} -SIGin{ii+1} );%GR_{i+1,i+1}
    prod_SIG_h = SIGi1_h{ii+1} * prod_SIG_h;
    new_GRi1 = Gii_next * prod_SIG_h;%GR_{i+1,1}
    %current inside the i-th slice
    Hij = sparse(triu(nearest00{2-mod(ii,2)}.id)) * (-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest00{2-mod(ii,2)}.xa - nearest00{2-mod(ii,2)}.xb).*(nearest00{2-mod(ii,2)}.ya + nearest00{2-mod(ii,2)}.yb)/2);
    current_intra(:,:,ii) = 2*imag( Hij .*(last_GRi1 * GAMMAlead_left * last_GRi1' ).' );
    
    current_inter(:,:,ii) = 2*imag( H01{-mod(ii,2)+2} .*(new_GRi1 * GAMMAlead_left * last_GRi1' ).' );
    last_GRi1 = new_GRi1;
end
ii = Nslice;
Hij = sparse(triu(nearest00{2-mod(ii,2)}.id)) * (-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest00{2-mod(ii,2)}.xa - nearest00{2-mod(ii,2)}.xb).*(nearest00{2-mod(ii,2)}.ya + nearest00{2-mod(ii,2)}.yb)/2);
current_intra(:,:,ii) = 2*imag( Hij .*(last_GRi1 * GAMMAlead_left * last_GRi1' ).' );
save local_current.mat current_intra current_inter
toc