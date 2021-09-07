function s = drawstates(y)
% Drawing from the 10 components indicators needed to construct the mixture distribution
% note that these are already demeaned!
% [prob mean variance]
% 10 componet mixture parameters Omori et al. 2004 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[prob			mean			variance]
pars10 = [ ...

0.04775    1.34744   0.17788;
0.13057    0.73504   0.26768;
0.20674    0.02266   0.40611;
0.22715   -0.85173   0.62699;

0.12047   -3.46788   1.57469;
0.05591   -5.55246   2.54498;
0.01575   -8.68384   4.16591;
0.00115  -14.65000   7.33342;];

pars  = pars10;
C     = size(pars,1);
seqc	= 1:C;  % sequnce of components
N 		= size(y,1);
% y   = randn(N,1);

% create probalities p(s=j|y*,h) prop. to p(y*|h,s)p(s=j)
p_		= normpdf(repmat(y,1,C),repmat(pars(:,2)',N,1),repmat(sqrt(pars(:,3))',N,1)).*repmat(pars(:,1)',N,1);
p			= p_./repmat(sum(p_,2),1,C);	% normalise so that probs sum to one.

% NOTE: carefull here, if normpdf is evaluated far from the mean, with probs
% close to zero, then p_ will be a vector of zeros and hence it will not be
% possible to draw from the mnrnd with p=0!

% 
msqc  = repmat(seqc,N,1)';   % create matrix of sequence numbers
sI    = mnrnd(1,p)==1;       % make indicator of states which are equal to 1.
s     = msqc(sI');           % states s=j with prob(s=j|y*,h)

