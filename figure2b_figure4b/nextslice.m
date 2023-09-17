function [nslice, update_res_sites] = nextslice(last_nu, res_sites)
%%% among the remaining lattice sites, find the nearest-hopping sites to
%%% last_nu
del1 = [0  1/sqrt(3)];
del2 = [1/2 -1/(2*sqrt(3))];
del3 = [-1/2 -1/(2*sqrt(3))];
temp_slice=repelem(last_nu(:,1:3),3,1);
temp_slice(:,1:2) = temp_slice(:,1:2) + kron(last_nu(:,3),[del1;del2;del3]);%all possible sites
temp_slice(:,3)=-temp_slice(:,3);%give labels to them, A or B
temp_nu = res_sites(ismembertol(res_sites(:,1:3),temp_slice,10^-6,'ByRows',true),:);
nslice = sortrows(temp_nu);
update_res_sites = res_sites(~ismembertol(res_sites,nslice,10^-9,'ByRows',true),:);
end