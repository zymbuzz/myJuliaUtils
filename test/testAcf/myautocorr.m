function acf=myautocorr(y,numLags)

y         = y - mean(y)
T         = length(y)

% function assumes that there are no missings

numLags=abs(numLags); % assuming that cov-stationary

acf = nan(numLags+1,1);

for j = 0:numLags
     
     cross   = y(1:end-j).*y(j+1:end);

     acf(j+1) = sum(cross)/T;
end

acf = acf./acf(1); % Normalize
end
