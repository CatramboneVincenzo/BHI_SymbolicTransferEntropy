function[vec] = Nquantization(x, n)
% take x and quantize it in n symbols according to min-max transformation
% IN THIS VERSION THE MIN AND MAX OF THE SIGNAL ITSELF IS CONSIDERED
if numel(x) ~= max(size(x))
    error('x must be a vector');
elseif size(x,1)<size(x,2)
    x = x';
end

x_comp = linspace(min(x),max(x),n+1);

vec = x >= x_comp;
vec = sum(vec,2);
vec(vec==(n+1)) = n;

end