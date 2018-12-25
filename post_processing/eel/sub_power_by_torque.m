
tt = linspace(0,2,num);
ss = linspace(0,1,nl_totl+1);
a_max = 11.41;
k = 2*pi/0.59; omega = 2*pi;
a = a_max.*exp(ss-1);

for nt = 1:num
    kappa(nt,:) = a.*sin(k*ss-omega*tt(nt));
    kappa_dot(nt,:) = -omega*a.*cos(k*ss-omega*tt(nt));
end

for nt = 1:num
    for nl = 1:nl_totl+1
%         dtheta(nt,nl) = sum(kappa_dot(nt,1:nl).*dss(1:nl));
%         power_fluid(nt,nl) = torque_fluid1(nt,nl)*dtheta(nt,nl);
%         power_inert(nt,nl) = torque_inert1(nt,nl)*dtheta(nt,nl);
%         power_tt(nt,nl) = Torque_totl(nt,nl)*dtheta(nt,nl);
        power_fluid(nt,nl) = torque_fluid1(nt,nl)*kappa_dot(nt,nl);
        power_inert(nt,nl) = torque_inert1(nt,nl)*kappa_dot(nt,nl);
        power_tt(nt,nl)    = Torque_totl(nt,nl)*kappa_dot(nt,nl);
      
    end
end
power_totl = power_fluid + power_inert;

% power_fluid_fine_f = func_smoothing_in_2D(power_fluid, num, nl_totl, refine1, refine2);
% power_inert_fine_f = func_smoothing_in_2D(power_inert, num, nl_totl, refine1, refine2);
% power_totl_fine_f  = func_smoothing_in_2D(power_tt, num, nl_totl, refine1, refine2);
% % 
low_v = -6e-3;
top_v = -low_v;
fn    = 'power by torque';
xs = linspace(0,1,nl_totl+1);
figure;
%imagesc(1:64,1:200,power_fluid_fine_f);
%imagesc(xs,xt,power_tt);
imagesc(xs,xt,power_tt);

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


power_positive = zeros(nt,nl_totl+1);
for nt = 1:num
    for nl = 1:nl_totl+1
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
for i=1:nl_totl+1
    Torque_e(:,i)=kk(i)*kappa(:,i);
end
Torque_sum=Torque_totl+Torque_e;
for nt = 1:num
    for nl = 1:nl_totl+1
        power_sum(nt,nl)    = Torque_sum(nt,nl)*kappa_dot(nt,nl); 
    end
end
power_sum_positive = zeros(nt,nl_totl+1);
for nt = 1:num
    for nl = 1:nl_totl+1
        if power_sum(nt,nl)<0
            power_sum_positive(nt,nl) = 0;
        else
            power_sum_positive(nt,nl) = power_sum(nt,nl);
        end
    end
end
Torque_totl_fine_fs = func_smoothing_in_2D(Torque_sum, num, nl_totl, refine1, refine2);
%power_totl_fine_fs  = func_smoothing_in_2D(power_sum, num, nl_totl, refine1, refine2);
low_v = -1.5e-4;
top_v = -low_v;
fn    = 'total torque';
xs = linspace(0,1,nl_totl+1);
xxx = linspace(0,1,198);
figure;
%imagesc(1:64,1:200,power_fluid_fine_f);
imagesc(xs,xt,Torque_sum);
%imagesc(xxs,xxt,power_sum);
hold on
% for ii = 30:6:165
% contour(xxx(ii-2:ii+2),xt,Torque_totl_fine_fs(:,ii-2:ii+2),[aa(ii),aa(ii)],'--k','LineWidth',3)
% hold on
% contour(xxx(ii-2:ii+2),xt,Torque_totl_fine_fs(:,ii-2:ii+2),[-aa(ii),-aa(ii)],'--k','LineWidth',3)
% hold on
% end
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
figure
        

ds=0.01;
for nt = 1:num
    power_h1(nt) = sum(power_tt(nt,:)*ds);
    power_h2(nt) = sum(power_positive(nt,:)*ds);
end
sum(power_h1*dt*10)/2
sum(power_h2*dt*10)/2