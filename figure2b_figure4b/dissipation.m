function [TeffRL, Voltage_real,energy_current, T_i, F_ii_ee] = dissipation()
tic
%=========%
[Sample, Stochastic] = Parameter();
%                |2|         |3|
%     ---===========---
%     -1-===========-4-
%     ---===========---
%                |6|         |5|
%%
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
%% surface Green's function
%=========%left and right contacts
[GAMMAlead_left,GAMMAlead_right, SIGMAlead_left, SIGMAlead_right]=lead(Sample.EF,Sample.phiz,Lattice);
%=========%
%=========% side leads
SIGMAsl = -1i/2 * Sample.GAMMAsl;
%=========%
siteA = Lattice.siteA;
siteB = Lattice.siteB;
sa = min(siteA(:,1));
edgeA = siteA(siteA(:,1)==sa,:);
edgeB = siteB(siteB(:,1)==sa,:);
uc = sortrows([edgeA;edgeB],2);
%% partition scheme
res_sites=[Lattice.siteA;Lattice.siteB];
N = size(res_sites,1);
nu={};
nu{1}=uc;
temp_N = size(nu{1},1);
N_nu=[];% number of sites in each slice
N_nu(1)=temp_N;
ii=1;
res_sites = res_sites(~ismembertol(res_sites,nu{1},10^-9,'ByRows',true),:);
while temp_N<N
    last_nu=nu{ii}; 
    [temp_nu, res_sites] = nextslice(last_nu, res_sites);
    nu{ii+1}=temp_nu;
    ii = ii+1;
    N_nu(ii) = size(temp_nu,1);
    temp_N=temp_N+N_nu(ii);
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
%% generate Hamiltonian
%=======%
H00=cell(Nslice,1);
H01 = cell(2,1);
nearest00 = cell(2, 1);
nearest00{1} = isnearest(nu{1},nu{1});
nearest01_1 = isnearest(nu{1}, nu{2});%(nu{1})by(nu{2})
H01{1} = sparse(nearest01_1.id*(-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest01_1.xa - nearest01_1.xb).*(nearest01_1.ya + nearest01_1.yb)/2));
%%%the second kind of slice
nearest00{2}  = isnearest(nu{2},nu{2});
nearest01_2 = isnearest(nu{2}, nu{3});%(nu{2})by(nu{3})
H01{2} = sparse(nearest01_2.id*(-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest01_2.xa - nearest01_2.xb).*(nearest01_2.ya + nearest01_2.yb)/2));
SIGMA_dissipation = zeros(Lattice.num_dissipation,1);
tem_num_diss = 0;
for ii = 1 : 1 : Nslice
    nii = N_nu(ii);
    tem_tem_diss = nnz(nu{ii}(:,5)==7);
    tem_SIGMA_dissipation = Stochastic.dissipation.str * ( nu{ii}(:,5)==7 );
    SIGMA_dissipation(tem_num_diss+1 : tem_num_diss + tem_tem_diss) = tem_SIGMA_dissipation(nu{ii}(:,5)==7);
    tem_num_diss = tem_num_diss + tem_tem_diss;
    H00{ii} = sparse( ...
        Sample.V0*eye(nii)...%on-site energy
        + Sample.gamma1*diag(nu{ii}(:,4))...%edge potential
        + Stochastic.disorder.str * spdiags(rand(nii,1) - 0.5, 0, nii, nii)...%aanderson disorder
        + Sample.SQUID.str * spdiags(exp(-vecnorm((nu{ii}(:,1:2) - Sample.SQUID.location),2,2)/Sample.SQUID.radius), 0, nii, nii)...%long-range Scattering barrier
        + (ii == 1) * SIGMAlead_left + (ii == Nslice) * SIGMAlead_right...%left and right contatcs
        + (ii ~= 1) * (ii ~= Nslice) * SIGMAsl * diag( nu{ii}(:,5)==2|nu{ii}(:,5)==3|nu{ii}(:,5)==5|nu{ii}(:,5)==6 )...%side leads
        + spdiags(tem_SIGMA_dissipation, 0, nii, nii)...%Buttiker probes
        + nearest00{2-mod(ii,2)}.id*(-Sample.t).*exp(1i*2*pi*Sample.phiz*(nearest00{2-mod(ii,2)}.xa - nearest00{2-mod(ii,2)}.xb).*(nearest00{2-mod(ii,2)}.ya + nearest00{2-mod(ii,2)}.yb)/2)...
        );
end
%% recursive Green's function
%%%initialize
gri_i = inv( Sample.EF*speye(N_nu(1)) - H00{1} );
grlead_i = gri_i;
gri_lead = gri_i;
grlead_lead = gri_i;
all_lead = nu{1}(:,[1,2,5,6]);%col1=x, col2=y, col3=lead, col4 = distribute_site
for ii = 2 : Nslice
    nii = N_nu(ii);

    tem_gri_i = gri_i;
    tem_grlead_i = grlead_i;
    tem_grlead_lead = grlead_lead;
    tem_gri_lead = gri_lead;
    tem_all_lead = all_lead;

    tem_gri_i = inv(Sample.EF*speye(nii) - H00{ii} - H01{mod(ii,2) + 1}' * tem_gri_i * H01{mod(ii,2) + 1});
    tem_grlead_i = tem_grlead_i * H01{mod(ii,2) + 1} * tem_gri_i;
    tem_grlead_lead = tem_grlead_lead + tem_grlead_i * H01{mod(ii,2) + 1}' * tem_gri_lead;
    tem_gri_lead = tem_gri_i * H01{mod(ii,2) + 1}' * tem_gri_lead;

    id_lead = (nu{ii}(:,5)>0)|(nu{ii}(:,6)>0);
    all_lead = [tem_all_lead;
        nu{ii}(id_lead,[1,2,5,6])];
    gri_i = tem_gri_i;
    grlead_lead = [tem_grlead_lead, tem_grlead_i(:,id_lead);
        tem_gri_lead(id_lead, :),tem_gri_i(id_lead, id_lead)];
    grlead_i = [tem_grlead_i;
        tem_gri_i(id_lead,:)];
    gri_lead = [tem_gri_lead, tem_gri_i(:, id_lead)];
end

tem2_all_lead = sortrows([all_lead, (1:size(all_lead,1)).'],[3,1,2]);%col1=x, col2=y, col3=lead, col4 = distribute_site, col5 = order for sort
GR = grlead_lead(tem2_all_lead(:,5), tem2_all_lead(:,5));
grlead_lead = GR( logical(tem2_all_lead(:,3)),logical(tem2_all_lead(:,3)) );

%% calculate transmission matrix
GAMMA_T = blkdiag( GAMMAlead_left, Sample.GAMMAsl * speye(2 * (Sample.sl_wid + 1)), GAMMAlead_right, Sample.GAMMAsl * speye(2 * (Sample.sl_wid + 1)), spdiags(2i*SIGMA_dissipation,0,Lattice.num_dissipation, Lattice.num_dissipation) );%for calculating transmission
SIGMA_T = blkdiag( SIGMAlead_left, SIGMAsl * speye(2 * (Sample.sl_wid + 1)), SIGMAlead_right, SIGMAsl * speye(2 * (Sample.sl_wid + 1)), spdiags(SIGMA_dissipation,0,Lattice.num_dissipation, Lattice.num_dissipation) );
temT = GAMMA_T*grlead_lead*GAMMA_T.*conj(grlead_lead);
sumM=blkdiag( ones(1, Sample.NWid), kron(speye(2),ones(1,Sample.sl_wid+1)), ones(1, Sample.NWid), kron(speye(2),ones(1,Sample.sl_wid+1)), speye(Lattice.num_dissipation) );
T = real(sumM*temT*sumM');
%% calculate Hall
V_R = 0;
V_L = 1;
Tphi = T([2,3,5,6, 7:end],[2,3,5,6, 7:end]);% except for the left and right leads
TRphi = T(4,[2,3,5,6, 7:end]);%here
TphiR = T([2,3,5,6, 7:end],4);%here
TphiL = T([2,3,5,6, 7:end],1);%here,
TLphi = T(1,[2,3,5,6, 7:end]);%here,
W = -Tphi + diag(sum(T([2,3,5,6, 7:end],:),2));
mu = V_R +  W\TphiL * (V_L - V_R);
Voltage = [V_L;mu(1:2);V_R;mu(3:4);mu(5:end)];
TeffRL = T(4,1) + TRphi*mu;
Voltage_real = Voltage(1:6);
%% calculate energy dissipation
[V_q, V_p] = meshgrid(Voltage,Voltage);
energy_current_p = sum(1/2 * T .* (V_p - V_q).^2, 2);% unit is e^2/h * (\delta V)^2
energy_current_p_tilde = [energy_current_p(1)/N_nu(1)*ones(N_nu(1),1);...
    energy_current_p(2)/(Sample.sl_wid+1)*ones((Sample.sl_wid+1),1);...
    energy_current_p(3)/(Sample.sl_wid+1)*ones((Sample.sl_wid+1),1);...
    energy_current_p(4)/N_nu(end)*ones(N_nu(end),1);...
    energy_current_p(5)/(Sample.sl_wid+1)*ones((Sample.sl_wid+1),1);...
    energy_current_p(6)/(Sample.sl_wid+1)*ones((Sample.sl_wid+1),1);...
    energy_current_p(7:end)];
%% calculate temperature 
T_bg = 0.1/1;%backgroud temperature,k_{B}T_{bg}/(e\delta V)
A = pi^2/6*(-T + diag(sum(T,2)));
id_dissipation_temper = false(Lattice.num_dissipation,1);
id_dissipation_temper(randperm(Lattice.num_dissipation, round(Stochastic.dissipation.temper* Lattice.num_dissipation))) = true;
id_temper = logical([zeros(6,1);id_dissipation_temper]);
A_tilde = A(id_temper, id_temper);
Ti_square = A_tilde\ (energy_current_p(id_temper)-T_bg^2*sum( A(id_temper, ~id_temper), 2 ));
Ti = zeros(Lattice.num_dissipation,1);
Ti(id_dissipation_temper) = (sqrt(Ti_square)-T_bg)/T_bg;
%%% record
all_sites = sortrows([[Lattice.siteA;Lattice.siteB],zeros(N,2)],[5,1,2]);
all_sites(end-length(energy_current_p_tilde)+1:end,7) = energy_current_p_tilde;
all_sites(end-Lattice.num_dissipation+1:end,8) = Ti;
all_sites = sortrows(all_sites);
energy_current = all_sites(:,7);
T_i = all_sites(:,8);
%% calculate energy distribution of electrons
distribut_site = tem2_all_lead(tem2_all_lead(:,4)>0, :);%col1=x, col2=y, col3=lead, col4 = distribute_site, col5 = order for sort
distribute_gr = full(GR(tem2_all_lead(:,4)>0, tem2_all_lead(:,3)>0));% distribute_site by all lead
distribute_site_gr = full(GR(tem2_all_lead(:,4)>0, tem2_all_lead(:,4)>0));% distribute_site by distribute_site
distribut_site  = sortrows([distribut_site(:,1:4), (1:size(distribut_site,1)).']);%col1=x, col2=y, col3=lead, col4 = distribute_site, col5 = order for sort
distribute_gr  = distribute_gr( logical(distribut_site(:,5)), :);
distribute_site_gr  = distribute_site_gr (logical(distribut_site(:,5)),logical(distribut_site(:,5)));
n_ii_ee = zeros(size(distribut_site , 1), length(Sample.distribute_EF));
ldos_ii_ee = zeros(size(distribut_site , 1), length(Sample.distribute_EF));
for ii = 1 : size(distribut_site , 1)
    ldos_ii_ee(ii,:) = -1/pi*imag( distribute_site_gr(ii, ii) );
    for ee = 1 : length(Sample.distribute_EF)
        SIG_in = blkdiag( (Sample.distribute_EF(ee)<Voltage(1))*GAMMAlead_left, (Sample.distribute_EF(ee)<Voltage(2))*Sample.GAMMAsl * speye((Sample.sl_wid + 1)), (Sample.distribute_EF(ee)<Voltage(3))*Sample.GAMMAsl * speye((Sample.sl_wid + 1)), (Sample.distribute_EF(ee)<Voltage(4))*GAMMAlead_right, (Sample.distribute_EF(ee)<Voltage(5))*Sample.GAMMAsl * speye((Sample.sl_wid + 1)), (Sample.distribute_EF(ee)<Voltage(6))*Sample.GAMMAsl * speye((Sample.sl_wid + 1)), spdiags((Sample.distribute_EF(ee)<Voltage(7:end)).*2i.*SIGMA_dissipation,0,Lattice.num_dissipation, Lattice.num_dissipation) );
        n_ii_ee(ii, ee) = 1/(2*pi) * distribute_gr(ii,:)*SIG_in*(distribute_gr(ii,:)');
    end
end
F_ii_ee = n_ii_ee./ldos_ii_ee;
toc
end