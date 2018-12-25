
for nt = 1:num
    for nl = 2:nl_body-1
        nslt = (nl-2)*nbody_peri+2:nl*nbody_peri+1;
        sts_int_s(nt,nl,1:3) = 0.5*sum(frc_body_t(nt,nslt,1:3))/ds(nl);
        sts_int_sl(nt,nl,1:3) = 0.5*sum(frc_body_t(nt,nslt,1:3));
    end
    nl = 1;
    sts_int_s(nt,1,1:3) = 2*(frc_body_t(nt,1,1:3)+0.5*sum(frc_body_t(nt,2:nbody_peri+1,1:3)))/(ds(nl)+ds(nl+1));
    sts_int_sl(nt,1,1:3) = frc_body_t(nt,1,1:3)+0.5*sum(frc_body_t(nt,2:nbody_peri+1,1:3));
    nl = nl_body;
    sts_int_s(nt,nl_body,1:3) = 2*(frc_body_t(nt,npoint_body,1:3)+0.5*sum(frc_body_t(nt,npoint_body-nbody_peri:npoint_body-1,1:3)))/(ds(nl-1)+ds(nl));
    sts_int_sl(nt,nl_body,1:3) = frc_body_t(nt,npoint_body,1:3)+0.5*sum(frc_body_t(nt,npoint_body-nbody_peri:npoint_body-1,1:3));
end

%smoothing the sts_int_s

numt = num*refine1;
xt = linspace(0,num, numt)/(num/2);
nums = (nl_totl-1)*refine2;
xs = linspace(0,nl_totl-1, nums)/(nl_totl-1);
for j = 1:2
    sts_temp = sts_int_s(:,:,j);
    sts_temp_f = func_smoothing_in_2D(sts_temp, num, nl_totl,refine1, refine2);
    sts_int_fine_f(1:numt,1:nums,j) = sts_temp_f;
end

%save force_fluid sts_int_s

for j = 1:2
     sts_temp = sts_int_s(:,:,j); 
     sts_temp_f = func_smoothing_in_2D(sts_temp, num, nl_totl, 1, 1);
     sts_int_sm(:,:,j) = sts_temp_f;
end


index = 2;
if index ==1
    low_v = -4e-3;
    top_v =  4e-3;
    fn = 'Fx';
elseif index ==2
    low_v = -5e-3;
    top_v =  5e-3;
    fn = 'Fy';
end

figure();
xxs = linspace(0,1,nl_totl);
xxt = linspace(0,2,num);
imagesc(xxs,xxt,sts_int_s(:,:,index));
hold on
contour(xs(3:end),xt,sts_int_fine_f(:,3:end,index),[1e-10,1e-10],'--k','LineWidth',3)
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution from CFD'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0.0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off
% hold on
% contour(xs(3:end),xt,sts_int_fine_f(:,3:end,index),[1e-80 1e-80],'linewidth',3,'linestyle','--','color','k')
