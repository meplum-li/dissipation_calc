function Lattice = hexagonal_lattice(Lx,Ly,pltshow)
%%% this function plot a graphene sized by Lx and Ly
%rectangle device
%Lx is a two-component vector [a, b], defining a interval of Lx
%Ly is a two-component vector [a, b], defining a interval of Ly
%we place the atom A at (0, 0), thus atom B locates at (0, 1/sqrt(3))
%bond length 1/sqrt(3)
[Sample, Stochastic] = Parameter();
edge_wid = Sample.edge_wid;
edge_wid_x = Sample.edge_wid_x;
startA = [0 0];% origin of first cell
dA = [0 0];% A relative to origin
dB = [0 1/sqrt(3)];% B relative to origin
%primitive vectors of a unit cell for both sublattice A and B
a1 = [1 0];
a2 = [1/2 sqrt(3)/2];
%% boundary for sublattice A
bdyA = ([a1' a2'])\([kron(Lx, [1 1]); kron([1 1], Ly)]-kron(ones(1,4),startA'));% boundary of sublattice A
[ARi, ARj] = meshgrid(floor(min(bdyA(1,:))):ceil(max(bdyA(1,:))), floor(min(bdyA(2,:))):ceil(max(bdyA(2,:))));
siteA=startA + ARi(:)*a1 + ARj(:)*a2;%generate A sites
idRA = (Lx(1)-startA(1)-siteA(:,1))<=10^-4& (siteA(:,1)-Lx(2)+startA(1))<=10^-4 & (Ly(1)-startA(2)-siteA(:,2))<=10^-4 & (siteA(:,2)-Ly(2)+startA(2))<=10^-4;
siteA = siteA(idRA,:);
%% for sublattice B
bdyB = ([a1' a2'])\([kron(Lx, [1 1]); kron([1 1], Ly)]-kron(ones(1,4),startA'+dB'));
[BRi, BRj] = meshgrid(floor(min(bdyB(1,:))):ceil(max(bdyB(1,:))), floor(min(bdyB(2,:))):ceil(max(bdyB(2,:))));
siteB=startA + dB + BRi(:)*a1 + BRj(:)*a2;
idRB = (Lx(1)-(startA(1))-siteB(:,1))<=10^-3 & (siteB(:,1)-Lx(2)+(startA(1)))<=10^-3 & (Ly(1)-(startA(2))-siteB(:,2))<=10^-3 & (siteB(:,2)-Ly(2)+(startA(2)))<=10^-3 ;
siteB = siteB(idRB,:);
%% bond between A and B
del1 = [0  1/sqrt(3)];
del2 = [1/2 -1/(2*sqrt(3))];
del3 = [-1/2 -1/(2*sqrt(3))];
tem_A = repelem(siteA,3,1);
A2B = tem_A + kron(ones(size(siteA,1),1),[del1;del2;del3]);
idA2B = ismembertol(A2B,siteB,10^-9,'ByRows',true);
bond = [tem_A(idA2B,:),A2B(idA2B,:)];%[Ax, Ay, Bx, By]
%% device show
Len = [min([siteA(:,1);siteB(:,1)]),max([siteA(:,1);siteB(:,1)])];
Wid = [min([siteA(:,2);siteB(:,2)]),max([siteA(:,2);siteB(:,2)])];
if pltshow
    figure('Visible',pltshow,'Position',[0 0 (Len(2)-Len(1))*50 (Wid(2)-Wid(1))*50])
    cc=rand(1,3);
    plot(siteA(:,1),siteA(:,2),'o','MarkerSize',4,'MarkerFaceColor',cc,'MarkerEdgeColor',cc,'DisplayName','Sublattice A');
    hold on
    cc=abs(1-cc);
    plot(siteB(:,1),siteB(:,2),'o','MarkerSize',4,'MarkerFaceColor',cc,'MarkerEdgeColor',cc,'DisplayName','Sublattice B');
    legend('AutoUpdate','off')
    for ii = 1 : size(bond,1)
        hold on
        plot(bond(ii,1:2:3),bond(ii,2:2:4),'-k','LineWidth',1)
    end
    xlim(Len)
    ylim(Wid)
    daspect([1 1 1])
end
%% other features of sites
Lattice.siteA = zeros(size(siteA, 1), 6);%col1=x,col2=y, col3=1(A)|-1(B), col4=1(edge_potential)|0, col5=lead_number(1,2,3,4,5,6 for real leads and 7 for fictious leads), col6=whether to calculate energy distribution
Lattice.siteB = zeros(size(siteB, 1), 6);
Lattice.siteA(:,1:3) = [siteA,ones(size(siteA,1),1)];
Lattice.siteB(:,1:3) = [siteB,-ones(size(siteB,1),1)];
real_Wid = Wid(2) - Wid(1);
real_Len = Len(2) - Len(1);
%======% for edge potential
Lattice.siteA(:,4) = ~( (Ly(1) + edge_wid*real_Wid <= Lattice.siteA(:,2))& (Lattice.siteA(:,2) <= Ly(2) - edge_wid*real_Wid)&(Lx(1) + edge_wid_x*real_Len <= Lattice.siteA(:,1)) );
Lattice.siteB(:,4) = ~( (Ly(1) + edge_wid*real_Wid <= Lattice.siteB(:,2))& (Lattice.siteB(:,2) <= Ly(2) - edge_wid*real_Wid)&(Lx(1) + edge_wid_x*real_Len <= Lattice.siteB(:,1)) );
%======%
%======%the fifth column to label contacts (1,2,3,4,5,6 for real leads and 7 for fictious leads)
%                |2|         |3|
%     ---===========---
%     -1-===========-4-
%     ---===========---
%                |6|         |5|
topsl_x = Sample.topsl_x;%x coordinate of the start point of each side lead at the top edge
topsl_y = Sample.topsl_y;
botsl_x = Sample.botsl_x;%x coordinate of the start point of each side lead at the bottom edge
botsl_y = Sample.botsl_y;
sl_wid = Sample.sl_wid;
Lattice.siteA(:,5) =Lattice.siteA(:,5) + (Lattice.siteA(:,1) == Lx(1));% contact 1
Lattice.siteB(:,5) =Lattice.siteB(:,5) + (Lattice.siteB(:,1) == Lx(1));% contact 1
Lattice.siteA(:,5) =Lattice.siteA(:,5) + 2 * ( (abs(Lattice.siteA(:,2) - Ly(2))<1e-6)&((Lattice.siteA(:,1)-topsl_x(1) - sl_wid) < 1e-6)&((topsl_x(1)-Lattice.siteA(:,1)) < 1e-6) );%  contact 2
Lattice.siteA(:,5) =Lattice.siteA(:,5) + 3 * ( (abs(Lattice.siteA(:,2) - Ly(2))<1e-6)&((Lattice.siteA(:,1)-topsl_x(2) - sl_wid) < 1e-6)&((topsl_x(2)-Lattice.siteA(:,1)) < 1e-6) );%  contact 3
Lattice.siteA(:,5) =Lattice.siteA(:,5) + 4 * (Lattice.siteA(:,1) == Lx(2));%  contact 4
Lattice.siteB(:,5) =Lattice.siteB(:,5) + 4 * (Lattice.siteB(:,1) == Lx(2));%  contact 4
Lattice.siteB(:,5) =Lattice.siteB(:,5) + 5 * ( (abs(Lattice.siteB(:,2) - Ly(1))<1e-6)&((Lattice.siteB(:,1)-topsl_x(2) - sl_wid) < 1e-6)&((topsl_x(2)-Lattice.siteB(:,1)) < 1e-6) );%  contact 5
Lattice.siteB(:,5) =Lattice.siteB(:,5) + 6 * ( (abs(Lattice.siteB(:,2) - Ly(1))<1e-6)&((Lattice.siteB(:,1)-topsl_x(1) - sl_wid) < 1e-6)&((topsl_x(1)-Lattice.siteB(:,1)) < 1e-6) );%  contact 6
%%% contatc 7, for Buttiker probes
A_possible_dissipation = find(Lattice.siteA(:,5)==0);
id_A_dissipation = randperm(length(A_possible_dissipation), round(length(A_possible_dissipation)*Stochastic.dissipation.dens));
Lattice.siteA(A_possible_dissipation(id_A_dissipation), 5) = 7;
B_possible_dissipation = find(Lattice.siteB(:,5)==0);
id_B_dissipation = randperm(length(B_possible_dissipation), round(length(B_possible_dissipation)*Stochastic.dissipation.dens));
Lattice.siteB(B_possible_dissipation(id_B_dissipation), 5) = 7;
% 
% Lattice.siteA([18051, 21859],5) = 7;
% 
Lattice.num_dissipation = nnz(Lattice.siteA(:,5)==7)+nnz(Lattice.siteB(:,5)==7);
%======%
%======%
%%% whether to calculate energy distribution
% find(Lattice.siteA(:,1)==180&abs(Lattice.siteA(:,2)-60.62)<0.1)
% Lattice.siteA([18129,37173,47973,  10961,29951, 45369, 3011, 18051,37119],6) = 1;
% Lattice.siteA([3011,5631,9051,13271,18051,22851,27651,32451,37119,41059],6) = 1;
%======%
Lattice.bond = bond;
Lattice.Len = Len;
Lattice.Wid = Wid;
fprintf("length of the graphene is %.2fa\n", Len(2)-Len(1))
fprintf("width of the graphene is %.2fa\n", Wid(2)-Wid(1))
end