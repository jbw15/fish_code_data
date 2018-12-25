load('dm.mat');
load('dm3_1.mat');
load('area.mat');
load('ss.mat');
%load test.mat

%smoothing the dm
% select = 1:3:64;
% dm_f = spline(select,dm3(select), 1:64);
% dm3=dm_f;

for nt = 1:num
    
    for nl = 1:nl_totl -2
        dx = hh(nt,nl+1,1) - hh(nt,nl,1);
        dh_dx(nt,nl) = (hh(nt,nl+1,2) - hh(nt,nl,2))/dx;
        dhdt_dx(nt,nl) = (vel_int_s(nt,nl+1,2) - vel_int_s(nt,nl,2))/dx;
    end

    nl = nl_totl-1;
    dhdt_dx(nt,nl) = 2*dhdt_dx(nt,nl-1) - dhdt_dx(nt,nl-2);
    dh_dx(nt,nl) = 2*dh_dx(nt,nl-1) - dh_dx(nt,nl-2);
    
    nl = nl_totl;
    dh_dx(nt,nl) = 2*dh_dx(nt,nl-1) - dh_dx(nt,nl-2);
    
    for nl = 2:nl_totl -2
        dx_2 = hhs(nt,nl+1,1) - hhs(nt,nl-1,1);       % ntbm
        d2hdx2(nt,nl)= (dh_dx(nt,nl+1)-dh_dx(nt,nl-1))/dx_2;
        dmdx(nt,nl) = (dm3(nl+1)-dm3(nl-1))/dx_2;
    end
    nl=1;
    d2hdx2(nt,nl) = 2*d2hdx2(nt,nl+1) -d2hdx2(nt,nl+2);
    dmdx(nt,nl) = 2*dmdx(nt,nl+1) -dmdx(nt,nl+2);
    
    nl=nl_totl-1;
    d2hdx2(nt,nl) = 2*d2hdx2(nt,nl-1) -d2hdx2(nt,nl-2);
    dmdx(nt,nl) = 2*dmdx(nt,nl-1) -dmdx(nt,nl-2);
end

for nt = 2:num-1
    accel_s(nt,:,1:2) = (vel_int_s(nt+1,:,1:2)-vel_int_s(nt-1,:,1:2))/2/dt/n_interval;
end
nt = 1;
accel_s(nt,:,1:2)  = 2*accel_s(nt+1,:,1:2) - accel_s(nt+2,:,1:2);
nt = num;
accel_s(nt,:,1:2)  = 2*accel_s(nt-1,:,1:2) - accel_s(nt-2,:,1:2);

% for nt = 1:num
%     dhdt_dx(nt,1)  =  2*dhdt_dx(nt,2)-dhdt_dx(nt,3);
%     dhdt_dx(nt,49) =  (dhdt_dx(nt,48)+dhdt_dx(nt,50))/2;
%     dhdt_dx(nt,59) =  2/3*dhdt_dx(nt,58)+1/3*dhdt_dx(nt,61); 
%     dhdt_dx(nt,60) =  1/3*dhdt_dx(nt,59)+2/3*dhdt_dx(nt,61); 
% end

temp = accel_s(:,:,2)+2*U*dhdt_dx+U^2*d2hdx2;
dm_t = repmat(dm3,num,1);

force_reactive1 = -dm_t.*temp;
force_reactive2 = -U*dmdx.*((vel_int_s(:,:,2)+U*dh_dx(:,1:nl_totl-1)));
force_reactive = force_reactive1 + force_reactive2;
% force_reactive = force_reactive1;
force_reactive_fine_f = func_smoothing_in_2D(force_reactive, num, nl_totl,refine1, refine2);

for nt=1:num  
    fy_rea(nt,:)=force_reactive(nt,:)./ds;
end
fy_rea_distri_fine_f = func_smoothing_in_2D(fy_rea, num, nl_totl, refine1, refine2);
rea_amp = sqrt(2*var(fy_rea_distri_fine_f));
xxs = linspace(0,1,nl_totl-1);
figure();
imagesc(xxs,xt, force_reactive);
axis xy;
colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
hold on
title('reactive force')
set(gca,'FontSize',20);
colormap('jet');
% caxis([-2e-2 2e-2])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')

% %break
% 
% figure()
% plot(xs,fy_rea_distri_fine_f(:,:)')
% hold on
% % p1 = plot(xs,res_amp,'linewidth',4,'linestyle','-.');
% % hold on
% % p2 = plot(xs,rea_amp,'linewidth',4,'linestyle',':');
% % hold on
% p3 = plot(xs,rea_amp,'r-','linewidth',4);
% title('amplitude')
% lgd = legend('Resistance','Reactive','CFD');
% set(lgd,'position',[0.18 0.75 0.24 0.15])
% xlabel('Head<-Position->Tail','FontName','Times','FontSize',20);
% ylabel('Lateral Force','FontName','Times','FontSize',20);
% axis tight
% %xlim([0,1])
% set(gca,'FontSize',20)
