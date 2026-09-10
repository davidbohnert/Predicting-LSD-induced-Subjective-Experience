function q = benjamini_hochberg(p)
%BENJAMINI_HOCHBERG False-discovery-rate adjusted p-values.
validateattributes(p, {'numeric'}, {'vector','real','finite','>=',0,'<=',1});
originalSize = size(p);
p = p(:);
[sortedP, order] = sort(p);
m = numel(p);
adjusted = sortedP .* m ./ (1:m)';
adjusted = flipud(cummin(flipud(adjusted)));
adjusted = min(adjusted, 1);
q = zeros(m, 1);
q(order) = adjusted;
q = reshape(q, originalSize);
end
