for nl = 2:nl_totl
    dss(nl) = ss(nl)-ss(nl-1);
end
dss(1) = 2*dss(2)-dss(3);

 for nt = 1:num
      power_hist1(nt) = sum(power_tt(nt,:).*dss);
%      power_hist1(nt) = sum(power_tt(nt,:));
 end
 
 for nt = 1:num
 power_totl1(nt,:) = power_tt(nt,:).*dss;
 end
 
for nt = 1:num
    for nl = 1:nl_totl-1
        power_fluid2(nt,nl)  =    sts_int_s(nt,nl,2)*vel_int_s(nt,nl,2)*ds(nl) ...
                                + sts_int_s(nt,nl,1)*vel_int_s(nt,nl,1)*ds(nl);
        power_inert2(nt,nl)  =   force_inert(nt,nl,2)*vel_int_s(nt,nl,2) ...
                                + force_inert(nt,nl,1)*vel_int_s(nt,nl,1);
                            
         power_fs2(nt,nl) =  power_fluid2(nt,nl)/ds(nl);
         power_is2(nt,nl) =  power_inert2(nt,nl)/ds(nl);
    end

end

power_totl2 = - power_fluid2 - power_inert2;
power_s2    = - power_fs2    - power_is2;

% ms_fine_f = linspace(0,1,126)
% % yyaxis left
ms = linspace(0,1,nl_totl-1);
plot(ss,sum(power_positive)/2,'b','LineWidth',4)
hold on
plot(ss,sum(power_tt)/2,'-.r','LineWidth',4)
hold on
plot(ss,sum(power_sum_positive)/2,':g','LineWidth',4)
hold on
plot(ms,sum(power_s2)/2,'--c','LineWidth',4)
set(gca,'ycolor','k');
title('Work')
ylabel('work','FontName','Times','FontSize',20);
% yyaxis right
% plot(ms,sum(-power_fluid2)/2,'-c','LineWidth',3)
axis tight
set(gca,'FontSize',20)
set(gca,'ycolor','k');
xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
% ylabel('work','FontName','Times','FontSize',20);
for nt = 1:num
     power_hist2(nt) = sum(power_totl2(nt,:));
%      power_hist2(nt) = sum(power_s2(nt,:));
end
sum(power_hist1*dt*10)/2
sum(power_hist2*dt*10)/2

