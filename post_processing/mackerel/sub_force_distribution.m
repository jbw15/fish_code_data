
for nt = 1:num
       for nl = 2:nl_body-1
           nslt = (nl-2)*nbody_peri+2:nl*nbody_peri+1;
           sts_int_s(nt,nl,1:3) = 0.5*sum(frc_body_t(nt,nslt,1:3))/ds(nl);
           sts_int_sl(nt,nl,1:3) = 0.5*sum(frc_body_t(nt,nslt,1:3));
       end
end
    
    nl =1;
    sts_int_s(nt,nl,1:3) = 2*(0.5*sum(frc_body_t(nt,2:nbody_peri+1,1:3)) ...
        + frc_body_t(nt,nl,1:3))/(ds(nl)+ds(nl+1));
    sts_int_sl(nt,nl,1:3) = 2*(0.5*sum(frc_body_t(nt,2:nbody_peri+1,1:3)) ...
        + frc_body_t(nt,nl,1:3));
    
    t_top = zeros(1,nl_tail-1);
    t_bot = zeros(1,nl_tail-1) + 100;
       
    for nt = 1:num
       
       for nn = 1:nl_tail-1
           kk          = 0;
           stress_temp = zeros(1,3);
           force_temp  = zeros(1,3);
           temp        = zeros(1,3);
           temp2       = zeros(1,3);
           
           for np = 1:nfine_tail
               px_tail =pnt_fine_tail(nt,np,1);
               if(px_tail>=hh(nt,nl_body+nn,1) && px_tail<hh(nt,nl_body+nn+1,1))
                   kk = kk+1;
                   temp(1:3)   = sts_fine_tail(nt,np,1:3);
                   temp2(1:3)  = frc_fine_tail(nt,np,1:3);
                   
                   stress_temp =  stress_temp + temp;
                   force_temp  =  force_temp  + temp2;
                   hz          =  abs(pnt_fine_tail(nt,np,3));
                   
                   if(hz > t_top(nn))
                       t_top(nn) = hz;
                   end
                   if(hz <t_bot(nn))
                       t_bot(nn) = hz;
                   end
               end
               
               ncount(nn) = kk;
               ht =(t_top(nn)-t_bot(nn))*2;
               if(nn==14)
                   ht = 0.005;
               end
               height(nn) = ht;
 
 %              sts_int_s(nt,nn+nl_body,1:3) = stress_temp(1:3)/kk*ht;
                sts_int_s(nt,nn+nl_body,1:3) = force_temp/ds(nn+nl_body);
                sts_int_sl(nt,nn+nl_body,1:3) = force_temp;
           end
       end
       sts_int_s(nt,nl_body,1:3) = (sts_int_s(nt,nl_body-1,1:3) ...
           + sts_int_s(nt,nl_body+1,1:3))/2;
       sts_int_sl(nt,nl_body,1:3) = (sts_int_sl(nt,nl_body-1,1:3) ...
           + sts_int_sl(nt,nl_body+1,1:3))/2;
    end

% smoothing the sts_int_s
% save force_fluid sts_int_sl
numt = num*refine1;
xt = linspace(0,num, numt)/(num/2);
nums = (nl_totl-1)*refine2;
xs = linspace(0,nl_totl-1, nums)/(nl_totl-1);
for j = 1:2
     sts_temp = sts_int_s(:,:,j); 
     sts_temp_f = func_smoothing_in_2D(sts_temp, num, nl_totl, refine1, refine2);
     sts_int_fine_f(1:numt,1:nums,j) = sts_temp_f;
end

for j = 1:2
     sts_temp = sts_int_s(:,:,j); 
     sts_temp_f = func_smoothing_in_2D(sts_temp, num, nl_totl, 1, 1);
     sts_int_sm(:,:,j) = sts_temp_f;
end

index = 2;
if index ==1
    low_v = -5e-3;
    top_v =  5e-3;
    fn = 'Fx';
elseif index ==2
    low_v = -2e-2;
    top_v =  2e-2;
    fn = 'Fy';
end
xxs = linspace(0,1,nl_totl-1);
xxt = linspace(0,2,num);
% the segment is 64 now
figure();
% imagesc(xs,xt,sts_int_fine_f(:,:,index));
imagesc(xxs,xxt,sts_int_s(:,:,index));
hold on

%contour(pot(3:end),xxt,sts_int_s(:,3:end,index),[0,0],'--k','LineWidth',3)
%imagesc(xxx,xt,sts_temp);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([fn  ' distribution from CFD'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')
% hold on
% contour(xs(3:end),xt,sts_int_s(:,3:end,index),[1e-5 1e-5],'linewidth',3,'linestyle','--','color','k')
