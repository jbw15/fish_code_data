% close all

dx = 0.01;
xx=0:dx:0.98;
nx=length(xx);
load('ss.mat')
ns=length(ss);
for i=1:length(ss)-1
    pot(i)=(ss(i+1)+ss(i))/2;
end
xx2=0.005:dx:0.99;
nx=length(xx2);
ns=length(pot);
for i=1:nx
    for j=1:ns
        if xx2(i)>pot(j)&&xx2(i)<=pot(j+1)
            ind2(i)=j;
            c2(i,1)=(pot(j+1)-xx2(i))/(pot(j+1)-pot(j));
            c2(i,2)=(xx2(i)-pot(j))/(pot(j+1)-pot(j));
            break
        end
    end
end
fy1 = fy_rea; fdy = fy_res; fmy = fy_cfd;
for i=1:nx
    k=ind2(i);
    FAY2(:,i)=c2(i,1)*fy1(:,k)+c2(i,2)*fy1(:,k+1);
    FDY2(:,i)=c2(i,1)*fdy(:,k)+c2(i,2)*fdy(:,k+1);
    FMY2(:,i)=c2(i,1)*fmy(:,k)+c2(i,2)*fmy (:,k+1);
end
A =FAY2;B =FDY2;C =FMY2;
% A = fy1; B = fdy; C = fmy;
a0=sqrt(2*var(A));b0=sqrt(2*var(B));c0=sqrt(2*var(C));
tt2=1:400;

for j=2:98
    rr1(j)=fminsearch(@(c)(sum((a0(j)*sin(2*pi/200*tt2+c)-1/3*(A(:,j-1)'+A(:,j)'+A(:,j+1)')).^2)),0.01);
    rr2(j)=fminsearch(@(c)(sum((b0(j)*sin(2*pi/200*tt2+c)-1/3*(B(:,j-1)'+B(:,j)'+B(:,j+1)')).^2)),0.01);
    rr3(j)=fminsearch(@(c)(sum((c0(j)*sin(2*pi/200*tt2+c)-1/3*(C(:,j-1)'+C(:,j)'+C(:,j+1)')).^2)),0.01);
end
rr1(1  )=fminsearch(@(c)(sum((a0(1  )*sin(2*pi/200*tt2+c)-A(:,1  )').^2)),0.01);
rr1(j+1)=fminsearch(@(c)(sum((a0(j+1)*sin(2*pi/200*tt2+c)-A(:,j+1)').^2)),0.01);
rr2(1  )=fminsearch(@(c)(sum((b0(1  )*sin(2*pi/200*tt2+c)-B(:,1  )').^2)),0.01);
rr2(j+1)=fminsearch(@(c)(sum((b0(j+1)*sin(2*pi/200*tt2+c)-B(:,j+1)').^2)),0.01);
rr3(1  )=fminsearch(@(c)(sum((c0(1  )*sin(2*pi/200*tt2+c)-C(:,1  )').^2)),0.01);
rr3(j+1)=fminsearch(@(c)(sum((c0(j+1)*sin(2*pi/200*tt2+c)-C(:,j+1)').^2)),0.01);
for j=2:98
    r(j,:)=fminsearch(@(c)(sum((c(1)*(A(:,j-1)+A(:,j)+A(:,j+1))+c(2)*(B(:,j-1)+B(:,j)+B(:,j+1))-sum(C(:,j-1:j+1),2)).^2)),[0.1 0.2]);
end
r(1  ,:)=fminsearch(@(c)(sum((c(1)*(A(:,1)+A(:,2))+c(2)*(B(:,1)+B(:,2))-sum(C(:,1:2),2)).^2)),[0.1 0.2]);
r(j+1,:)=fminsearch(@(c)(sum((c(1)*(A(:,end-1)+A(:,end))+c(2)*(B(:,j)+B(:,j+1))-sum(C(:,j:j+1),2)).^2)),[0.1 0.2]);

for j=86:93
    r(j,1)=fminsearch(@(cc1)(sum((cc1*(A(:,j-1)+A(:,j)+A(:,j+1))-sum(C(:,j-1:j+1),2)).^2)),0.1);
    r(j,2)=fminsearch(@(cc2)(sum((cc2*(B(:,j-1)+B(:,j)+B(:,j+1))-sum(C(:,j-1:j+1),2)).^2)),0.1);
end
for j=1:99
    A3(:,j)=r(j,1)*A(:,j);
    B3(:,j)=r(j,2)*B(:,j);
end
C3=A3+B3;

figure
plot([0 1],[1 1],'Color',[0.8 0.8 0.8],'linewidth',3)
% hold on
% plot([0 1],[0 0],'Color',[0.8 0.8 0.8],'linewidth',3)
hold on
xx1 = [0.855 0.915 0.915 0.855];yy1 = [-2,-2,10,10];
H_F1 = fill(xx1,yy1,[0.8 0.8 0.8]);
set(H_F1,{'LineStyle'},{'none'}) 
plot(0.005:dx:0.845,r(1:85,1),'linewidth',4,'linestyle','-.','color','b')
hold on
plot(0.855:dx:0.915,r(86:92,1),'linewidth',4,'linestyle','-.','color','b')
hold on
plot(0.925:dx:0.985,r(93:99,1),'linewidth',4,'linestyle','-.','color','b')
hold on
plot(0.005:dx:0.845,r(1:85,2),'linewidth',4,'linestyle',':','color','r');
hold on
plot(0.855:dx:0.915,r(86:92,2),'linewidth',4,'linestyle',':','color','r');
hold on
plot(0.925:dx:0.985,r(93:99,2),'linewidth',4,'linestyle',':','color','r');
title('weight')
%set(h,'Box','off')
% plot([0.855 0.855],[min(r(:,2)),max(r(:,2))])
% plot([0.935 0.935],[min(r(:,2)),max(r(:,2))])
% % legend('c1','c2')
AA=sqrt(2*var(A3));BB=sqrt(2*var(B3));CC=sqrt(2*var(C));
xlabel('Head<-Body->Tail');
axis([0 1 0 10])
xlim([0,1])
set(gca,'FontSize',20)
figure
xx3 = [0.855 0.915 0.915 0.855];yy3 = [min(sqrt(2*var(B3(:,1:85)))),min(sqrt(2*var(B3(:,1:85)))),max(sqrt(2*var(A3(:,86:99)))),max(sqrt(2*var(A3(:,86:99))))];
H_F1 = fill(xx3,yy3,[0.8 0.8 0.8]);
set(H_F1,{'LineStyle'},{'none'}) 
hold on
plot(0.005:dx:0.845,AA(1:85),'b','linestyle','-.','linewidth',4)
hold on
plot(0.855:dx:0.915,AA(86:92),'b','linestyle','-.','linewidth',4)
hold on
plot(0.925:dx:0.985,AA(93:99),'b','linestyle','-.','linewidth',4)
hold on
plot(0.005:dx:0.845,BB(1:85),'r','linestyle',':','linewidth',4)
hold on
plot(0.855:dx:0.915,BB(86:92),'r','linestyle',':','linewidth',4)
hold on
plot(0.925:dx:0.985,BB(93:99),'r','linestyle',':','linewidth',4)
hold on
plot(0.005:dx:0.845,CC(1:85),'k','linewidth',4)
hold on
plot(0.855:dx:0.915,CC(86:92),'k','linewidth',4)
hold on
plot(0.925:dx:0.985,CC(93:99),'k','linewidth',4)
xlabel('Head<-Position->Tail','FontName','Times','FontSize',25);
axis tight
xlim([0,1])
ylim
title('amplitude')
% xticks([0 1])
set(gca,'FontSize',20)



rr12=rr1;rr12(1:41)=rr1(1:41)+2*pi;%you can change the range according to your need



rr22=rr2;rr22(1:27)=rr2(1:27)+2*pi;
rr32=rr3;rr32(1:38)=rr3(1:38)+2*pi;
% rr32 = rr3; rr32(1:23) = rr3(1:23)+2*pi;
figure
xx4 = [0.855 0.915 0.915 0.855];yy4 = [2*pi,2*pi,-2*pi,-2*pi];
H_F1 = fill(xx4,yy4,[0.8 0.8 0.8]);
set(H_F1,{'LineStyle'},{'none'}) %������ɫ���߿�
hold on
plot(xx2,rr12,'linewidth',4,'color','b','linestyle','-.')
hold on
plot(xx2,rr22,'linewidth',4,'color','r','linestyle',':')
hold on
plot(xx2,rr32,'linewidth',4,'linestyle','-','color','k')
% legend('EBT','Drag','CFD')
title('phase')
xlim([0,1])
xticks([0 1])
ylim([-2*pi 2*pi])
yticks([-2*pi -pi 0 pi 2*pi])
yticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
% xlabel('Head<-Position->Tail','FontName','Times','FontSize',25)
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
% ylabel('\fontsize{20}\fontname{Times new roman}work')
set(gca,'FontSize',20);

set(gca,'ydir','reverse')
FMY2_temp_f = func_smoothing_in_2D(FMY2, num, 100, refine1, refine2);
FMY2_int_fine_f(1:numt,1:198) = FMY2_temp_f;
xxx = linspace(0.005,0.99,198);
%============change the segment from 64 to 100 to get the figure of Fy====%
%=======in order to compare with phase figure=============================% 
figure();
low_v=-0.02;
top_v = -low_v;
%imagesc(xs,xt,sts_int_fine_f(:,:,index));
imagesc(xx2,xxt,FMY2);
% hold on
% contour(xxx,xt,FMY2_int_fine_f,[0,0],'--k','LineWidth',3)
%imagesc(xxx,xt,sts_temp);
axis xy;
cb = colorbar;
xlabel('Head<-Position->Tail','FontName','Times');
ylabel('Time','FontName','Times');
title([' Fy from CFD'], 'Fontname', 'Times', 'FontSize', 20)
set(gca,'FontSize',20,'xtick',[0.0 0.5 1.0], 'ytick',[0 1.0 2.0]);
colormap('jet');
set(cb,'xtick',[low_v 0 top_v])
caxis([low_v top_v])
xlabel('\fontsize{20}\fontname{Times new roman}Head<-Position->Tail')
ylabel('\fontsize{20}\fontname{Times new roman}Time')