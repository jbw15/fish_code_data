
ds = 0.01;

for nt = 1:num
    power_hist1(nt) = sum(power_tt(nt,:).*ds);
end

 for nt = 1:num
 power_totl1(nt,:) = power_tt(nt,:).*ds;
 end
 
 for nt = 1:num
     for nl = 1:nl_totl
%         power_fluid2(nt,nl)  =    sts_int_s(nt,nl,2)*vel_int_s(nt,nl,2)*ds(nl) ...
%                                 + sts_int_s(nt,nl,1)*vel_int_s(nt,nl,1)*ds(nl);
         power_fluid2(nt,nl) = sts_int_sl(nt,nl,2)*vel_int_s(nt,nl,2) ...
                             + sts_int_sl(nt,nl,1)*vel_int_s(nt,nl,1);
        power_inert2(nt,nl)  =   force_inert(nt,nl,2)*vel_int_s(nt,nl,2) ...
                                + force_inert(nt,nl,1)*vel_int_s(nt,nl,1);
                            
        power_fs2(nt,nl) =  power_fluid2(nt,nl)/ds;
        power_is2(nt,nl) =  power_inert2(nt,nl)/ds;
    end

 end
power_totl2 = - power_fluid2 - power_inert2 ; 
power_s2 = -power_fs2 - power_is2;
ms = linspace(0,1,nl_totl);
% yyaxis left
plot(ss,sum(power_positive)/2,'b','LineWidth',4)
hold on
plot(ss,sum(power_tt)/2,'-.r','LineWidth',4)
hold on
plot(ss,sum(power_sum_positive)/2,':g','LineWidth',4)
hold on
plot(ms,sum(power_s2)/2,'c--','LineWidth',4)
title('work')
%legend('positive work','total work','total work(muscle+elasticity)','work done to fluid')
legend('W^{+}','W','W^{+}_{elasticity}','W_{fluid}')
set(gca,'ycolor','k');
ylabel('work','FontName','Times','FontSize',20);
% yyaxis right
% plot(ms,sum(-power_fs2)/2,'-c','LineWidth',3)
axis tight
set(gca,'FontSize',20)
set(gca,'ycolor','k');
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);


for nt = 1:num
    power_hist2(nt) = sum(power_totl2(nt,:));
end
sum(power_hist1*dt*10)/2
sum(power_hist2*dt*10)/2
% %==========================================================================
% 
% for nt = 1:num
%     power_totl3(nt,:) = vel_body_t(nt,:,1).*frc_body_t(nt,:,1) + vel_body_t(nt,:,2).*frc_body_t(nt,:,2);
% end
% for nt = 1:num
%     power_totl31(nt) = sum(power_totl3(nt,:));
%     power_totl32(nt) = sum(power_inert2(nt,:));
%     power_hist3(nt) = -power_totl31(nt) - power_totl32(nt);
% end
% 
% % for nt = 1:num
% %     power_hist3(nt) = sum(power_totl3_1(nt,:));
% % end
% xp = (1:400)/200;
% figure()
% axes('FontName','Times','FontSize',18)
% hold on
% h1 = plot(xp, power_hist1,'k-','LineWidth',2)
% hold on
% h2 = plot(xp, power_hist2,'r-','LineWidth',3)
% % hold on
% % h3 = plot(xp, power_hist3,'b-','LineWidth',3)
% xlim([0 2])
% % ylim([-1e-4 7e-4])
% 
% xlabel('t/T', 'FontName','Times','FontSize',18)
% ylabel('total power', 'FontName','Times','FontSize',18)
% legend([h1 h2], 'By torque', 'By force')
% % legend([h1 h2], 'By torque', 'By force1')
% hold off
% 
% 
low_v = -2e-5;
top_v = -low_v;
fn    = 'power by force';
xs = linspace(0,1,100);
figure;
imagesc(xs,xt,power_totl2);
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
