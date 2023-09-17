function [Sample, Stochastic] = Parameter()
%=======%Fermi energy
% Sample.EF=linspace(-2.28,2.28,2);
Sample.EF = 0.05;
%=======%
Sample.t = 2.7;
Sample.V0=0;
%=======%geometry setting of device
Sample.gamma2 = 0;
Sample.M = 0;% stagger potential
Sample.EFlead=1.5*Sample.t;
Sample.eta = 10^(-6);
Sample.phiz = -0.004;%gauge choice: A=(-By,0,0), phiz = eB(a_0)^2/h, e<0
Sample.B_real = -Sample.phiz*2*2.0678*10^(-15)/(0.246*10^(-9))^2;
Sample.Lx = [0,202+1/2];% length in the center of startA along x
Sample.Ly = [-sqrt(3)/2*1/3, sqrt(3)* (120 - 0.5) ];% length in the center of startA along y。here, Ly = [-sqrt(3)/2*(n1+1/3),sqrt(3)/2*n2]，n1=n2+1
%=======%
%edge potential
Sample.gamma1 = 0.1;%strength
Sample.edge_wid = 60/(Sample.Ly(2) - Sample.Ly(1));%ratio of the edge width to the device width in y direction
Sample.edge_wid_x = 20/(Sample.Lx(2) - Sample.Lx(1));%width in x direction
%=======%
Sample.NWid=2*round((Sample.Ly(2) + sqrt(3)/2)/sqrt(3));%Width of device in y direction
Sample.NLen = ((Sample.Lx(2)-Sample.Lx(1)))/(0.5)+1;%zigzag, length of device in x direction
Sample.N = Sample.NLen*Sample.NWid;%total number of lattice sites
%%% geometry setting of side contacts
%                |2|         |3|
%     ---===========---
%     -1-===========-4-
%     ---===========---
%                |6|         |5|
Sample.topsl_x = [30, 170] + 0.5;%x coordinate of the start point of each side lead at the top edge
Sample.topsl_y = Sample.Ly(2);
Sample.botsl_x = Sample.topsl_x;%x coordinate of the start point of each side lead at the bottom edge
Sample.botsl_y = Sample.Ly(1);
Sample.sl_wid = 9;% width of each side contact
Sample.GAMMAsl  = 2*pi*2;
%%% SQUID setting
% potential function: Sample.SQUID.str * exp(-|r - r_0|/Sample.SQUID.radius)
Sample.SQUID.str = -0.13*Sample.t;
Sample.SQUID.radius = 10;%radius is about 2.5nm
% Sample.SQUID. Wid= ( Sample.Ly(2) - Sample.Ly(1) ) * Sample.edge_wid*0.8;
% Sample.SQUID.Len = 10;
Sample.SQUID.location = [mean(Sample.topsl_x + [Sample.sl_wid, 0]), Sample.Ly(1) + ( Sample.Ly(2) - Sample.Ly(1) ) * (Sample.edge_wid/2)];
%%% setting for calculating energy distribution
Sample.distribute_EF = linspace(-0.1,1.1,100);
%============%setting for dissipation and disorder
Stochastic.dissipation.str = -1i/2*0.1*Sample.t;%strength of dissipation%%SIGMA
Stochastic.dissipation.dens = 0.005;%density of dephasing sites
Stochastic.dissipation.temper = 0.1;%for Buttiker probes，the ratio of probes which increase temperature
Stochastic.disorder.str = 0*Sample.t;%strength of disorder, [-disorder.str/2,disorder.str/2]
Stochastic.disorder.dens = 1;%density of disorder
Stochastic.aver = 4*320;%number of different realizations of randomness
Stochastic.core =32;% number of core for your computer
end