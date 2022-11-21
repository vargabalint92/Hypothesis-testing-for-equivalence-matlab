function [is_accept,p_crit,p_act] = sign_rank_equivalence(d, qpl1,qpl2,alpha)
%% Implementation of the paired-data analogue of the Mann-Whitney-Wilcoxon
% test for equivalence of continuous distributions.
% The continuity assumption relates to the intraindividual
% differences D_i. For details
% see Wellek S (2010) Testing statistical hypotheses of equivalence and
% noninferiority. Second edition,Par. 5.4.
% d - row vector with the intraindividual differences for all n pairs as components



qplct = (qpl1+qpl2)/2;
eps =(qpl2-qpl1)/2;

u = 0;
n = length(d);

for i=1:(n-1)
    for j=(i+1):n
        u = u + floor(0.5*(sign(d(i)+d(j))+1));
    end
end
zeta = 0;

for i =1:(n-2)
    for j=(i+1):(n-1)
        for k=(j+1):n
            zeta = zeta + floor(0.5*(sign(min(d(i)+d(j),d(i)+d(k))) + 1)) + ...
                floor(0.5*(sign(min(d(j)+d(i),d(j)+d(k))) + 1)) + ...
                floor(0.5*(sign(min(d(k)+d(i),d(k)+d(j))) + 1));
        end
    end
end


u = u*2/n/(n-1);
zeta = zeta*2/n/(n-1)/(n-2) - u^2;
sigmah = sqrt( (4*(n-2)*zeta + 2*u*(1-u) ) /n/(n-1) );
p_crit = sqrt(ncx2inv(alpha,1,(eps/sigmah)^2));

p_act = abs((u-qplct)/sigmah);

if p_act >= p_crit
    is_accept = 0;
else
    is_accept = 1;
end
end