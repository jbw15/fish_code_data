
for nt = 1:num
    
    %re-capture the center linew
    for nl = 2:nl_body
        nslt = (nl-2)*nbody_peri+2:(nl-1)*nbody_peri+1;
        hh(nt,nl,1:3) = mean(pnt_body_t(nt,nslt,1:3));
       
    end
     hh(nt,1,1:3)=pnt_body_t(nt,1,1:3);
     hh(nt,nl_body+1,1:3)=pnt_body_t(nt,npoint_body,1:3);

    
    for nl = 2:nl_body-1
        nslt = (nl-2)*nbody_peri+2:nl*nbody_peri+1;
        hhs(nt,nl,1:3) = mean(pnt_body_t(nt,nslt,1:3));
    end
    hhs(nt,1,1:3) = mean(pnt_body_t(nt,1:nbody_peri+1,1:3));
    hhs(nt,nl_body,1:3) = mean(pnt_body_t(nt,npoint_body-nbody_peri:npoint_body,1:3));
    %calculate ds of each element
    ss=zeros(1,nl_body);
    for nl = 1:nl_body
        dsx = hh(nt,nl+1,1) - hh(nt,nl,1);
        dsy = hh(nt,nl+1,2) - hh(nt,nl,2);
        ds(nl) = sqrt(dsx^2 + dsy^2);
        % calculate the center of each element
       % hhs(nt,nl,1:2) = (hh(nt,nl+1,1:2) + hh(nt,nl,1:2))/2;
    end
    %ds(1) = 2*ds(2) - ds(3);
    ss(1) = ds(1);
    for nl = 2:nl_body
        ss(nl) = ss(nl-1) + ds(nl);
    end
end

% save center_position hhs
