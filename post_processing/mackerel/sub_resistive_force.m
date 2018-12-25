
%============calculate the resistive force=================
temp_a = zeros(1,3);
temp_v = zeros(1,3);
for nt=1:num
    for nl=1:nl_totl-1
        temp_a(1:2) = accel_s(nt,nl,1:2);
        temp_v(1:2) = vel_int_s(nt,nl,1:2);
        temp_body(1:2) = hh(nt,nl+1,1:2) - hh(nt,nl,1:2);
        temp_a(3) = 0;
        temp_v(3) = 0;
        temp_body(3) = 0;
        accel_unit = temp_a(1:2)/norm(temp_a);
        veloc_unit = temp_v(1:2)/norm(temp_v);
        acc_amp = norm(temp_a);
        body_unit = temp_body(1:2)/norm(temp_body);
        if dot(veloc_unit,[-body_unit(2),body_unit(1)])>0
            vy_unit = [-body_unit(2),body_unit(1)];
        else 
            vy_unit = [body_unit(2),-body_unit(1)];
        end
        cos_theta = dot(veloc_unit,body_unit)/(norm(veloc_unit)*norm(body_unit));
        f_inert(nt,nl)  = acc_amp*(-dm(nl));               %inertial force
        force_inert(nt,nl,:)=f_inert(nt,nl)*accel_unit;
        
        if nl<=nl_body
            Cd = 1;
        else
            Cd = 1.0;
        end
        
%         fd(nt,nl)=-1/2*Cd*area(nl)*sum(temp_v.^2);          %resistive force
%         force_resist2 (nt,nl,:)=fd(nt,nl)*veloc_unit;
        fd(nt,nl)=-1/2*Cd*area(nl)*sum(temp_v.^2)*(1-cos_theta^2);          %resistive force
        force_resist2 (nt,nl,:)=fd(nt,nl)*vy_unit;
    end
    
end

% lateral force

Fy_resist2 =  force_resist2(:,:,2);
% Fy_cfd     =  sts_int_sm(:,:,2);
% Fx_cfd     =  sts_int_sm(:,:,1);
Fy_cfd     =  sts_int_s(:,:,2);
Fx_cfd     =  sts_int_s(:,:,1);
%  density of lateral force
for nt=1:num  
    fy_rea(nt,:)=force_reactive(nt,:)./ds;
    fy_res(nt,:)=Fy_resist2(nt,:)./ds;
    fy_cfd(nt,:)=Fy_cfd(nt,:);
end

fy_res_distri_fine_f = func_smoothing_in_2D(fy_res, num, nl_totl, refine1, refine2);
fy_rea_distri_fine_f = func_smoothing_in_2D(fy_rea, num, nl_totl, refine1, refine2);
fy_cfd_distri_fine_f = func_smoothing_in_2D(fy_cfd, num, nl_totl, refine1, refine2);
fy_EBT= fy_res_distri_fine_f + fy_rea_distri_fine_f;


low_v = -0.01;
top_v = -low_v;
fn  = 'fy ebt';
figure();
imagesc(xs,xt*2,fy_res_distri_fine_f);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 0.5 1.0 1.5 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
hold off


% hold on
% contour(xs(3:end),xt,sts_int_s(:,3:end,index),[1e-5 1e-5],'linewidth',3,'linestyle','--','color','k')

