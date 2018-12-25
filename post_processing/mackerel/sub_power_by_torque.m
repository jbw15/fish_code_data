
scale = 1;
tt = linspace(0,2,num);

a0 = 1*scale; a1= -3.2*scale; a2 = 5.6*scale;
k = 2*pi/1.0; omega = 2*pi;
A0 = a0 + a1*ss + a2*ss.^2;

for nt = 1:num
        kappa(nt,:)    = A0.*sin(k*ss - omega*tt(nt));
        kappa_dot(nt,:)= -omega*A0.*cos(k*ss - omega*tt(nt));
end

for nt = 1:num
    for nl = 1:nl_totl
%         power_fluid(nt,nl)  =  torque_fluid1(nt,nl)*kappa_dot(nt,nl);
%         power_inert(nt,nl)  =  torque_inert1(nt,nl)*kappa_dot(nt,nl);
        power_tt(nt,nl)     =  Torque_totl(nt,nl)*kappa_dot(nt,nl);
    end
end

%power_totl = power_fluid + power_inert;
% 
%  power_fluid_fine_f = func_smoothing_in_2D(power_fluid, num, nl_totl, refine1, refine2);
%  power_inert_fine_f = func_smoothing_in_2D(power_inert, num, nl_totl, refine1, refine2);
 power_totl_fine_f  = func_smoothing_in_2D(power_tt, num, nl_totl, refine1, refine2);
% 
low_v = -2e-3;
top_v = -low_v;
fn    = 'power';

figure;
%imagesc(1:64,1:200,power_fluid_fine_f);
imagesc(xss,xxt,power_tt);
hold on
contour(xs(33:113),xt,Torque_totl_fine_f(:,33:113),[1e-50 1e-50],'--k','LineWidth',3)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off

power_positive = zeros(nt,nl_totl);
for nt = 1:num
    for nl = 1:nl_totl
        if power_tt(nt,nl)<0
            power_positive(nt,nl) = 0;
        else
            power_positive(nt,nl) = power_tt(nt,nl);
        end
    end
end
T0 = sqrt(2*var(Torque_totl));
kapa0=sqrt(2*var(kappa));
kk=0.4*T0./kapa0;
for i=1:nl_totl
    Torque_e(:,i)=kk(i)*kappa(:,i);
end
Torque_sum=Torque_totl+Torque_e;
for nt = 1:num
    for nl = 1:nl_totl
        power_sum(nt,nl)    = Torque_sum(nt,nl)*kappa_dot(nt,nl); 
    end
end
power_sum_positive = zeros(nt,nl_totl);
for nt = 1:num
    for nl = 1:nl_totl
        if power_sum(nt,nl)<0
            power_sum_positive(nt,nl) = 0;
        else
            power_sum_positive(nt,nl) = power_sum(nt,nl);
        end
    end
end
Torque_totl_fine_fs = func_smoothing_in_2D(Torque_sum, num, nl_totl, refine1, refine2);
% power_totl_fine_fs  = func_smoothing_in_2D(power_sum, num, nl_totl, refine1, refine2);
low_v = -3e-4;
top_v = -low_v;
fn    = 'total torque';

figure;
%imagesc(1:64,1:200,power_fluid_fine_f);
%imagesc(xs,xt,power_tt);
%imagesc(xss,xxt,power_sum);
imagesc(xss,xxt,Torque_sum);
% imagesc(xs,xt,Torque_totl_fine_fs );
% hold on
% for ii = 30:6:95
% contour(xs(ii-2:ii+2),xt,Torque_totl_fine_fs(:,ii-2:ii+2),[aa(ii),aa(ii)],'--k','LineWidth',3)
% hold on
% contour(xs(ii-2:ii+2),xt,Torque_totl_fine_fs(:,ii-2:ii+2),[-aa(ii),-aa(ii)],'--k','LineWidth',3)
% hold on
% end
hold on
contour(xs(27:101),xt,Torque_totl_fine_fs(:,27:101),[0,0],'--k','LineWidth',3)
% hold on
% contour(xs(27:101),xt,Torque_totl_fine_fs(:,27:101),[-6.6150e-5,-6.6150e-5],'--k','LineWidth',3)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
%title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off
 
for nl = 2:nl_totl
    dss(nl) = ss(nl)-ss(nl-1);
end
dss(1) = 2*dss(2)-dss(3);
for nt = 1:num
    power_h1(nt) = sum(power_tt(nt,:).*dss);
    power_h2(nt) = sum(power_positive(nt,:).*dss);
end
sum(power_h1*10*dt)/2
sum(power_h2*10*dt)/2
% mean(power_h1)
% mean(power_h2)
