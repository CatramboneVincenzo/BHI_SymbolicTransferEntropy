function[vec] = tuples_conversion(x, l)
% take symbols vector x  of length N and convert it in a N-2 long vector of tuples of length l
vec = zeros(length(x)-l+1, l);
for ii = 1:l
    vec(:,ii) = x(ii:length(x)-l+ii);
end


