function [y] = moving_average(x, n)
y = nan(size(x));
% for i = (n+1):size(x, 2)-n
%     y(:, i) = mean(x(:, i-n:i+n), 2, 'omitnan');
% end
for i = 1:size(x, 2)
    y(:, i) = mean(x(:, max(i-n, 1):min(i+n, size(x,2))), 2, 'omitnan');
end
end