function [Sample, Stochastic] = Parameter()
%=======%费米能
% Sample.EF=linspace(-2.28,2.28,2);
Sample.EF = 0.05;
%=======%不需要修改
Sample.t = 2.7;
Sample.V0=0;
%=======%装置的几何设置
Sample.gamma2 = 0;
Sample.M = 0;% stagger potential
Sample.EFlead=1.5*Sample.t;
Sample.eta = 10^(-6);
Sample.phiz = -0.004;%规范 A=(-By,0,0), phiz = eB(a_0)^2/h, e<0
Sample.B_real = -Sample.phiz*2*2.0678*10^(-15)/(0.246*10^(-9))^2;
Sample.Lx = [0,202+1/2];% length in the center of startA along x
Sample.Ly = [-sqrt(3)/2*1/3, sqrt(3)* (120 - 0.5) ];% length in the center of startA along y。这里Ly = [-sqrt(3)/2*(n1+1/3),sqrt(3)/2*n2]，n1=n2+1
%=======%
%边缘势场
Sample.gamma1 = 0.1;%边缘势场的强度
Sample.edge_wid = 60/(Sample.Ly(2) - Sample.Ly(1));%上下边缘势场存在的宽度，即单侧外场宽度占整个Ly宽度的比例。取值范围：(0,0.5)；绝对长度大致保持在60。
Sample.edge_wid_x = 20/(Sample.Lx(2) - Sample.Lx(1));%左边边缘势垒存在的宽度，即外场宽度占整个Lx宽度的比例。取值范围：(0,0.5)；绝对长度大致保持在40。
%=======%
Sample.NWid=2*round((Sample.Ly(2) + sqrt(3)/2)/sqrt(3));%最左的格点数目，其中A和B子格子的格点交错排列；
Sample.NLen = ((Sample.Lx(2)-Sample.Lx(1)))/(0.5)+1;%zigzag构型，最下边一列（y=0&&y=-1/(2sqrt(3))）的A+B格点的数目；
Sample.N = Sample.NLen*Sample.NWid;%total number of lattice sites
%%% 侧面电极的几何设置
%                |2|         |3|
%     ---===========---
%     -1-===========-4-
%     ---===========---
%                |6|         |5|
Sample.topsl_x = [30, 170] + 0.5;%x coordinate of the start point of each side lead at the top edge
Sample.topsl_y = Sample.Ly(2);
Sample.botsl_x = Sample.topsl_x;%x coordinate of the start point of each side lead at the bottom edge
Sample.botsl_y = Sample.Ly(1);
Sample.sl_wid = 9;% 侧面电极的宽度，x方向格点数目-1
Sample.GAMMAsl  = 0*2*pi*2;
%%% SQUID设置
% potential function: Sample.SQUID.str * exp(-|r - r_0|/Sample.SQUID.radius)
Sample.SQUID.str = -0.13*Sample.t;
Sample.SQUID.radius = 10;%半径大约是2.5nm
% Sample.SQUID. Wid= ( Sample.Ly(2) - Sample.Ly(1) ) * Sample.edge_wid*0.8;
% Sample.SQUID.Len = 10;
Sample.SQUID.location = [mean(Sample.topsl_x + [Sample.sl_wid, 0]), Sample.Ly(1) + ( Sample.Ly(2) - Sample.Ly(1) ) * (Sample.edge_wid/2)];
%%% 测量电子分布的设置
Sample.distribute_EF = linspace(-0.1,1.1,100);
%============%耗散和无序的设置
Stochastic.dissipation.str = -1i/2*0.1*Sample.t;%strength of dissipation%%SIGMA
Stochastic.dissipation.dens = 0;%density of dephasing sites
Stochastic.dissipation.temper = 0.4;%虚拟导线中，不导热（升温）的电极比例
Stochastic.disorder.str = 0*Sample.t;%strength of disorder, [-disorder.str/2,disorder.str/2]
Stochastic.disorder.dens = 1;%density of disorder
Stochastic.aver = 1;%number of different realizations of randomness
Stochastic.core =32;%服务器的核心数
end