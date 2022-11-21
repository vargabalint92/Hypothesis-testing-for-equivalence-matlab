function [is_accept, p_crit, p_act] = mann_whitney_equivalence(x,y,eps1_,eps2_,alpha)
%% Equivalence Testing
%Implementation of the asymptotically distribution-free test for 
% equivalence of two continuous distributions in terms of the 
% Mann-Whitney-Wilcoxon functional. For details see Wellek S (2010) 
% Testing statistical hypotheses of equivalence and noninferiority. 
% Second edition, Par. 6.2.


m = length(x);
n = length(y);
eqctr = 1 + (eps2_-eps1_)/2 ;
eqleng = eps1_ + eps2_;

wxy = 0;
pihxxy = 0;
pihxyy = 0;


for i=1:m
    for j=1:n
        wxy = wxy + floor(0.5*(sign(x(i) - y(j)) + 1));
    end
end


for i=1:m
    for j1=1:(n-1)
        for j2=(j1+1):n
            pihxyy = pihxyy + floor(0.5*(sign(x(i) - max(y(j1),y(j2))) + 1));
        end
    end
end
       
for i1=1:(m-1)
    for i2=(i1+1):m
        for j=1:n
            pihxxy = pihxxy + floor(0.5*(sign(min(x(i1),x(i2)) - y(j)) + 1));
        end
    end
end


wxy = wxy/(m*n);
pihxxy  = pihxxy*2 / (m*(m-1)*n);
pihxyy  = pihxyy*2 / (n*(n-1)*m);
sigmah  = sqrt((wxy-(m+n-1)*wxy^2+(m-1)*pihxxy+(n-1)*pihxyy)/(m*n));

p_crit = sqrt(ncx2inv(alpha,1,(eqleng/2/sigmah)^2));
%p_crit = sqrt(pdf('Normal',x,mu,sigma));
p_act = abs((wxy-eqctr)/sigmah);
if (p_act >= p_crit) 
    is_accept = 1;
elseif (p_act < p_crit)  
    is_accept = 0;
else
    is_accept = 0;
end
