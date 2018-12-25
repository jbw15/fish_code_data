function [ ff_fine_f ] = func_smoothing_in_2D(ff, num, nl_totl, refine1, refine2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% smoothing the sts_int_s
width = 5;
for nl = 1:nl_totl-1
        for nt = 1:num
            if nt>=width+1 && nt <=num-width
                ff_smt(nt,nl) = mean(ff(nt-width:nt+width,nl));
            elseif nt < width +1
                ff_smt(nt,nl) = mean(ff(1:nt,nl));
            elseif nt> num-width
                ff_smt(nt,nl) = mean(ff(nt:num,nl));
            end
        end
end
width = 1;
for nt = 1:num
        for nl = 1:nl_totl-1
            if nl>=width+1 && nl <=nl_totl-1-width
                ff_sm(nt,nl) = mean(ff_smt(nt,nl-width:nl+width));
            elseif nl < width +1
                ff_sm(nt,nl) = mean(ff_smt(nt,1:nl));
            elseif nl > nl_totl-1-width
                ff_sm(nt,nl) = mean(ff_smt(nt,nl:nl_totl-1));
            end
        end
end

numt = num*refine1;
xt = linspace(1,num, numt);
nums = (nl_totl-1)*refine2;
xs = linspace(1,nl_totl-1, nums);

ff_fine =  zeros(numt,nl_totl-1);
ff_fine_f =  zeros(numt,nums);

for nl = 1:nl_totl-1
        temp2 = spline(1:num,ff_sm(:,nl), xt) ;
        ff_fine(1:numt,nl) = temp2;
end

for nt = 1:numt
      temp2 = spline(1:nl_totl-1, ff_fine(nt,:), xs) ;
      ff_fine_f(nt,1:nums) = temp2;
end

end

