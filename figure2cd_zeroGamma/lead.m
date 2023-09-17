function [GAMMAlead_left,GAMMAlead_right, SIGMAlead_left, SIGMAlead_right]=lead(EF,phiz,Lattice)
%%% 左右电极的自能
[Sample, Stochastic] = Parameter();
t = Sample.t;
V0=Sample.V0;
M = Sample.M;
gammaM = 0;
gamma1 = Sample.gamma1;
edge_wid = Sample.edge_wid;
gamma2 = Sample.gamma2;
% SIGMAlead = zeros(size(llead,1)+size(rlead,1));
eta = Sample.eta;
%% left
% unit cell for lead
siteA = Lattice.siteA;
siteB = Lattice.siteB;
sa_x = min(siteA(:,1));%A格点的横坐标最小值
edgeA = siteA(siteA(:,1)==sa_x|siteA(:,1)==sa_x+0.5,:);%最左边的元胞内的A格点
edgeB = siteB(siteB(:,1)==sa_x|siteB(:,1)==sa_x+0.5,:);%最左边的元胞内的B格点
uc = sortrows([edgeA;edgeB],2);%按照格点的y坐标，从小到大对最左边的元胞内的格点进行排序。uc的第三个分量标志子格子的类型，1 for A and -1 for B %unit cell
puc = uc - [1,zeros(1, size(uc,2) - 1)];%previous unit cell: first slice of left lead
nearest01 = isnearest(puc, uc);
LH01 = nearest01.id*(-t).*exp(1i*2*pi*phiz*(nearest01.xa - nearest01.xb).*(nearest01.ya + nearest01.yb)/2);
nearest00 = isnearest(puc, puc);
LH00 = V0*eye(size(puc,1))+M*diag(puc(:,3)).*exp(-gammaM*min([max(puc(:,2))-puc(:,2),puc(:,2)- min(puc(:,2))],[],2 ))...
    +gamma1*diag(uc(:,4))...
    +nearest00.id*(-t).*exp(1i*2*pi*phiz*(nearest00.xa - nearest00.xb).*(nearest00.ya + nearest00.yb)/2);%生成哈密顿量
%%% 开始迭代计算
ai = LH01';
bi = LH01;
ei = LH00;
eg = LH00;
zz = (EF+ 1i*eta)*speye(size(puc,1));
ss = zeros(50,1);
for ind = 1 : 50
    mm = inv(zz - ei);
    eg = eg + ai * mm * bi;
    ei = ei + ai * mm * bi + bi * mm * ai;
    ai = ai * mm * ai;
    bi = bi * mm * bi;
    ss(ind) = sum(sum(abs(ai) + abs(bi)));
    if ss(ind)<1e-14
        break
    end
end
gs_left = inv(zz - eg);%surface green function of lead left
edgeA = siteA(siteA(:,1)==sa_x,:);%最左边的A格点
edgeB = siteB(siteB(:,1)==sa_x,:);%最左边的B格点
llead = sortrows([edgeA;edgeB],2);%按照格点的y坐标，从小到大对最左边的格点进行排序。uc的第三个分量标志子格子的类型，1 for A and -1 for B %unit cell
nearest_lc = isnearest(puc,llead);%(uc)by(llead), llead指中心区最左边跟左电极接触的格点
H_lc = nearest_lc.id*(-t).*exp(1i*2*pi*phiz*(nearest_lc.xa - nearest_lc.xb).*(nearest_lc.ya + nearest_lc.yb)/2);
SIGMAlead_left = H_lc'*gs_left* H_lc;
%% right lead
% unit cell for right lead
la_x = max(siteA(:,1));%A格点的横坐标最大值
edgeA = siteA(siteA(:,1)==la_x|siteA(:,1)==la_x-0.5,:);%最右边的元胞内的A格点
edgeB = siteB(siteB(:,1)==la_x|siteB(:,1)==la_x-0.5,:);%最右边的元胞内的B格点
uc = sortrows([edgeA;edgeB],2);%1 for A and -1 for B %unit cell
nuc = uc + [1,zeros(1, size(uc,2) - 1)];%first slice of left lead
nearest01 = isnearest(nuc, uc);%(nuc)by(uc)
RH01 = nearest01.id*(-t).*exp(1i*2*pi*phiz*(nearest01.xa - nearest01.xb).*(nearest01.ya + nearest01.yb)/2);%(nuc)by(uc)
nearest00 = isnearest(nuc, nuc);
RH00 = V0*eye(size(nuc,1))...
    +gamma1*diag(nuc(:,4))...
    +nearest00.id*(-t).*exp(1i*2*pi*phiz*(nearest00.xa - nearest00.xb).*(nearest00.ya + nearest00.yb)/2);%生成哈密顿量
ai = RH01';
bi = RH01;
ei = RH00;
eg = RH00;
zz = (EF+ 1i*eta)*speye(size(nuc,1));
ss = zeros(50,1);
for ind = 1 : 50
    mm = inv(zz - ei);
    eg = eg + ai * mm * bi;
    ei = ei + ai * mm * bi + bi * mm * ai;
    ai = ai * mm * ai;
    bi = bi * mm * bi;
    ss(ind) = sum(sum(abs(ai) + abs(bi)));
    if ss(ind)<1e-14
        break
    end
end
gs_right = inv(zz - eg);%side lead right
edgeA = siteA(siteA(:,1)==la_x,:);
edgeB = siteB(siteB(:,1)==la_x,:);
rlead = sortrows([edgeA;edgeB],2);%按照格点的y坐标，从小到大对最右边的格点进行排序。uc的第三个分量标志子格子的类型，1 for A and -1 for B %unit cell
nearest_rc = isnearest(nuc,rlead);%(nuc)by(rlead), rlead指中心区最右边跟右侧电极接触的格点
H_rc = nearest_rc.id*(-t).*exp(1i*2*pi*phiz*(nearest_rc.xa - nearest_rc.xb).*(nearest_rc.ya + nearest_rc.yb)/2);% right lead to central region
SIGMAlead_right = H_rc'*gs_right* H_rc;%假设lead和central region的耦合跟导线内部一样
%%
GAMMAlead_left=1i*(SIGMAlead_left-SIGMAlead_left');
GAMMAlead_right = 1i*(SIGMAlead_right-SIGMAlead_right');
end