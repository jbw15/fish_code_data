load('dm.mat')
for nt = 1:num
    for nl = 2:nl_body-1
        nslt = (nl-2)*nbody_peri+2:nl*nbody_peri+1;
        hhs(nt,nl,1:3) = mean(pnt_body(nt,nslt,1:3));
    end
    hhs(nt,1,1:3) = mean(pnt_body(nt,1:nbody_peri+1,1:3));
    hhs(nt,nl_body,1:3) = mean(pnt_body(nt,npoint_body-nbody_peri:npoint_body,1:3));
    cx(nt)=sum(dm.*hhs(nt,:,1))/sum(dm);
    cy(nt)=sum(dm.*hhs(nt,:,2))/sum(dm);
end
(cx(400)-cx(1))/(dt*400*10)
(cy(400)-cy(1))/(dt*400*10)
