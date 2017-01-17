function torque = CalculateTorque(rho, DIGELocs, pow);
% Takes in rho, the correlation matrix for a network to all genes
% And DIGElocs, the position of the DIGE targets
% Calculate the torque of the difference in CDFs
% First, calculate center of mass, c
% Then, multiply c*A
corr_vec = rho(:);
[vals pos] = sort(corr_vec); % default: least to greatest

samp_vec = rho(:, DIGELocs);
samp_vec = samp_vec(:);

[a b c] = intersect(vals, samp_vec);

missed = ones(length(vals),1);
missed(b) = 0; % Arrange by order of 'vals'
hits = zeros(length(vals),1);
hits(b) = 1; %abs(vals(b)).^pow;

Pmiss = cumsum(missed./sum(missed));
Phit = cumsum(hits./sum(hits));
diff_cdf = Phit - Pmiss;

cut = min(find(vals>=0)); % Evaluate area on each side of rho=0
F.neg = diff_cdf(1:cut);
F.pos = diff_cdf(cut+1:end);

c_of_m.neg = sum(F.neg.*vals(1:cut))/sum(F.neg);
c_of_m.pos = sum(F.pos.*vals(cut+1:end))/sum(F.pos);

% Weight = area = Riemann sum
weight.neg = sum(F.neg.*[diff(vals(1:cut+1))]);
weight.pos = sum(F.pos.*[diff(vals(cut:end))]);
torque = (c_of_m.neg*weight.neg + c_of_m.pos*weight.pos); % Should be negative value and negative value, for greatest bimodality

