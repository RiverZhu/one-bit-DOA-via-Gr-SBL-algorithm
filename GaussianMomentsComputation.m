function [mout, vout] = GaussianMomentsComputation(y, tauin, phatin, vpin, vnin)
% GaussianMomentsComputation returns posterior mean and variance for E(z|y)
% and Var(z|y), where y=sign(z+tau+w)
% Input:
% - y: sign measurements
% - phatin: prior mean of z
% - vpin:  prior variance of z
% - vnin: additive noise of w
% - tauin: thresholds
% Output:
% - mout: E(Z | Y = y)
% - vout: Var(Z | Y = y)

% E(Z | Y = y)
alpha = (phatin+tauin) ./ sqrt(vpin+vnin);


C = y.*alpha;
CDF = normcdf(C,0,1);
ll = log(CDF);
                    %Find bad values that cause log cdf to go to infinity
I = find(C < -30);
                    %This expression is equivalent to log(normpdf(C)) for
                    %all negative arguments greater than -38 and is
                    %numerically robust for all values smaller.  DO NOT USE
                    %FOR LARGE positive x.
ll(I) = -log(2)-0.5*C(I).^2+log(erfcx(-C(I)/sqrt(2)));
temp = exp(log(normpdf(alpha))-ll)./sqrt(vpin+vnin);
% temp = normpdf(0,phatin+tauin,sqrt(vpin+vnin))./ normcdf(y.*alpha);
% temp = normpdf(alpha)./ (0.5*erfc(-y.*alpha/sqrt(2)))./sqrt(vpin+vnin);
mout = phatin+sign(y).*vpin.*temp;
% Var(Z | Y = y)
% temp1 = (phatin./vpin-tauin./vnin)./(1./vpin+1./vnin)-phatin;
temp1 = (phatin.*vnin-tauin.*vpin)./(vpin+vnin)-phatin;
if(isnan(temp1))
    display('warning, temp1 is nan');
end
vout = vpin+sign(y).*temp.*temp1.*vpin-(vpin.*temp).^2;

% max(vout)

% % E(Z | Y = y)
% alpha = (phatin+tauin) ./ sqrt(vpin+vnin);
% temp = normpdf(0,phatin+tauin,sqrt(vpin+vnin))./ normcdf(y.*alpha);
% mout = phatin+sign(y).*vpin.*temp;
% % Var(Z | Y = y)
% temp1 = (phatin./vpin-tauin./vnin)./(1./vpin+1./vnin)-phatin;
% vout = vpin+sign(y).*temp.*temp1.*vpin-(vpin.*temp).^2;
end

