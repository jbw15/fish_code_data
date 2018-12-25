load('dm.mat')    
for nt = 1:num
        
        % re-capture the center line
        for nl=2:nl_body
            nslt = (nl-2)*nbody_peri+2:(nl-1)*nbody_peri+1;
            hh(nt,nl,1:2) = mean(pnt_body(nt,nslt, 1:2));
        end
            hh(nt,1,1:2) = pnt_body(nt,1, 1:2);

        % tail before forking
        for nl=nl_body+1:nl_body+nl_tail_1-3
            nslt = (nl-nl_body+1)*ntail_vtcl+8;
            hh(nt,nl,1:2) = pnt_tail(nt,nslt, 1:2);
        end
            nl = nl_body+nl_tail_1-2;
            nslt = (nl-nl_body+1)*ntail_vtcl+9;
            hh(nt,nl,1:2) = pnt_tail(nt,nslt, 1:2);
        % tail after forking
         for nl = nl_body+nl_tail_1-1:nl_totl
            hh(nt,nl,1:2) =  pnt_tail(nt,nl+100, 1:2);
         end
           
      % calculate ds of each element
      for nl = 1:nl_totl-1
         dsx = hh(nt,nl+1,1) - hh(nt,nl,1);
         dsy = hh(nt,nl+1,2) - hh(nt,nl,2);
         ds(nl) = sqrt(dsx^2+dsy^2); 
         hhs(nt,nl,1:2) = (hh(nt,nl+1,1:2) + hh(nt,nl,1:2))/2;
      end
         ds(1) = 2*ds(2) - ds(3);
             cx(nt)=sum(dm.*hhs(nt,:,1))/sum(dm);
    cy(nt)=sum(dm.*hhs(nt,:,2))/sum(dm);
end

(cx(400)-cx(1))/(dt*400*10)
(cy(400)-cy(1))/(dt*400*10)