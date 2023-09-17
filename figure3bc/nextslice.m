function [nslice, update_res_sites] = nextslice(last_nu, res_sites)
%%% 在剩余未分配的格点中，寻找last_nu的最近邻的格点，作为下一个slice
del1 = [0  1/sqrt(3)];
del2 = [1/2 -1/(2*sqrt(3))];
del3 = [-1/2 -1/(2*sqrt(3))];
temp_slice=repelem(last_nu(:,1:3),3,1);
temp_slice(:,1:2) = temp_slice(:,1:2) + kron(last_nu(:,3),[del1;del2;del3]);%所有可能的最近邻节点
temp_slice(:,3)=-temp_slice(:,3);%所有可能的最近邻节点，标记A or B分类。这里，互为最近邻的格点，他们的格点类型是相反的。
temp_nu = res_sites(ismembertol(res_sites(:,1:3),temp_slice,10^-6,'ByRows',true),:);
nslice = sortrows(temp_nu);%默认从第一列开始，从小到大；如果第一列相同，则比较第二列
update_res_sites = res_sites(~ismembertol(res_sites,nslice,10^-9,'ByRows',true),:);
end