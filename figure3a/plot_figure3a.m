clearvars
load('local_current.mat')
[Sample, Stochastic] = Parameter();
Lattice = hexagonal_lattice(Sample.Lx,Sample.Ly,0);
flow=sortrows( [mean(Lattice.bond(:,[1,3]),2), mean(Lattice.bond(:,[2,4]),2),Lattice.bond] );
flow_direct = flow(:,[5,6])-flow(:,[3,4]);
sign_flow_direct = sign(flow_direct);
sign_flow_direct = (sign_flow_direct(:,1)==0).*sign_flow_direct(:,2) + (sign_flow_direct(:,1)~=0).*sign_flow_direct(:,1);
flow_direct = flow_direct .* sign_flow_direct;
flow_direct = flow_direct(flow_direct(:,1)~=0,:);
unit_flow_direct_x = flow_direct(:,1)./vecnorm(flow_direct,2,2);
unit_flow_direct_y = flow_direct(:,2)./vecnorm(flow_direct,2,2);
unit_flow_direct_x = reshape(unit_flow_direct_x,240,[]);
unit_flow_direct_y = reshape(unit_flow_direct_y,240,[]);
%%%for the final column
unit_flow_direct_x = [unit_flow_direct_x, zeros(240,1)];
unit_flow_direct_y = [unit_flow_direct_y, zeros(240,1)];

current_intra = squeeze(sum(current_intra,2));
current_inter = -squeeze(sum(current_inter,2));
current_x = zeros(240, 406);
current_y = zeros(240, 406);
for ii = 1 : 406
    if ii == 1
        current_x(:, ii) = current_inter(:, ii).*unit_flow_direct_x(:,ii);
        current_y(:, ii) = current_inter(:, ii).*unit_flow_direct_y(:,ii) + current_intra(:, ii);
    elseif ii == 406
            current_x(:, ii) = current_inter(:, ii-1).*unit_flow_direct_x(:, ii-1);
            current_y(:, ii) = current_inter(:, ii-1).*unit_flow_direct_y (:, ii-1)+ current_intra(:,ii);
    else
        current_x(:, ii) = current_inter(:, ii).*unit_flow_direct_x(:,ii) + current_inter(:, ii-1).*unit_flow_direct_x(:, ii-1);
        current_y(:, ii) = (current_inter(:, ii).*unit_flow_direct_y(:,ii) + current_inter(:, ii-1).*unit_flow_direct_y (:, ii-1))*3;
    end
end
all_sites = sortrows([Lattice.siteA;Lattice.siteB]);
loc_x = reshape(all_sites(:,1),[],Sample.NLen);
loc_y = reshape(all_sites(:,2),[],Sample.NLen);
% average
ncurrent_x = zeros(239, 203);
ncurrent_y = zeros(239, 203);
nloc_x = zeros(239, 203);
nloc_y = zeros(239, 203);
 for ii = 1 : 239
     if mod(ii, 2)==1
        ncurrent_x(ii,:) = [1 1]* current_x(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), -1), [1;0]));
        ncurrent_y(ii,:) = [1 1]* current_y(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), -1), [1;0])); 
        nloc_x(ii,:) = ([1 1]* loc_x(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), -1), [1;0])))./[ones(1,202)*6,4];
        nloc_y(ii,:) = ([1 1]* loc_y(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), -1), [1;0])))./[ones(1,202)*6,4];
     else
        ncurrent_x(ii,:) = [1 1]* current_x(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), 1), [0;1]));
        ncurrent_y(ii,:) = [1 1]* current_y(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), 1), [0;1]));
        nloc_x(ii,:) = ([1 1]* loc_x(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), 1), [0;1])))./[4,ones(1,202)*6];
        nloc_y(ii,:) = ([1 1]* loc_y(ii:(ii+1),:) *   (kron(eye(203), [1;1]) + kron(diag(ones(1,202), 1), [0;1])))./[4,ones(1,202)*6];
     end
 end
% current = -cat(1,squeeze(sum(current_intra,2)),squeeze(sum(current_inter,2)));
% current_y = current(1:240,:) + current(241:end,:) .*unit_flow_direct_y;
% current_x = current(241:end,:) .*unit_flow_direct_x;
angle = (ncurrent_x>0).*atand(ncurrent_y./ncurrent_x) + (ncurrent_x<0).*(180+atand(ncurrent_y./ncurrent_x)) + (ncurrent_x==0).*((ncurrent_y>0)*90 + (ncurrent_y<0)*270);
amplitude = sqrt(ncurrent_x.^2+ncurrent_y.^2);
% amplitude = abs(ncurrent_x) + abs(ncurrent_y);
x_ratio = 1:8:203;
y_ratio = 1:12:239;
figure
h=quiver(nloc_x(y_ratio,x_ratio), nloc_y(y_ratio,x_ratio),ncurrent_x(y_ratio,x_ratio),ncurrent_y(y_ratio,x_ratio),5,'LineWidth',2,'Color','w');
axis equal
xlim(Sample.Lx+[-2,2])
ylim(Sample.Ly+[-2,2])
hold on
surf(nloc_x,nloc_y,amplitude/10.2);%normalized the current
% colorbar
shading interp;
colormap jet;
view(0,90);
c=colorbar;
c.Label.String = 'J';
clim([0,0.1])
set(gca, 'FontSize', 20);
print('local_electric_current','-dpng','-r200')