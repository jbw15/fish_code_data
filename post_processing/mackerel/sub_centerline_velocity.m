for nt = 2:num-1
     if nt>1 && nt < num
        for j=1:2
            vel_int_s(nt,:,j) = (hhs(nt+1,:,j)-hhs(nt-1,:,j))/(2*dt*n_interval);
        end
     end
 end

 for j = 1:2
     vel_int_s(1,:,j) = 2*vel_int_s(2,:,j) - vel_int_s(3,:,j);
     vel_int_s(num,:,j) = 2*vel_int_s(num-1,:,j) - vel_int_s(num-2,:,j);
 end

%sub_original_curvature;
%vel_diff = vel_int_s(:,:,2) - kappa0_dot1;

figure();
imagesc(xs,xt,vel_int_s(:,:,2));
axis xy;
colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
hold on
% plot(Data001(:,1),Data001(:,2),'--k','linewidth',3)
% plot(Data002(:,1),Data002(:,2),'--k','linewidth',3)
% yticks([0.0025 1.0025  2.0025])
% yticklabels({0 1 2})
% xticks([-0.005 1.005])
% xticklabels({0 1})
title('lateral velocity distribution')
set(gca,'FontSize',20);
colormap('jet');
% caxis([-2e-1 2e-1])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
