function [th_max,h] = perfprof(A,th_max,colors,strings,linewidth)
%PERFPROF Performance profile.
% [th_max, h] = PERFPROF(A,th_max) produces a
% peformance profile for the data in the M-by-N matrix A,
% where A(i,j) > 0 measures the performance of the j’th solver
% on the i’th problem, with smaller values of A(i,j) denoting
% "better". For each solver theta is plotted against the
% probability that the solver is within a factor theta of
% the best solver over all problems, for theta on the interval
% [1, th_max].
% Set A(i,j) = NaN if solver j failed to solve problem i.
% TH_MAX defaults to the smallest value of theta for which
% all probabilities are 1 (modulo any NaN entries of A).
% h is a vector of handles to the lines with h(j)
% corresponding to the j’th solver.
    minA = min(A,[],2);
    if nargin < 2, th_max = max( max(A,[],2)./minA ); end
    tol = sqrt(eps); % Tolerance.
    [m,n] = size(A);
    for j = 1:n
        % m problems, n solvers.
        % Loop over solvers.
        col = A(:,j)./minA;
        % Performance ratios.
        col = col(~isnan(col)); % Remove NaNs.
        if isempty(col), continue; end
        theta = unique(col)';
        % Unique elements, in increasing order.
        r = length(theta);
        prob = sum( col(:,ones(r,1)) <= theta(ones(length(col),1),:) ) / m;
        % Assemble data points for stairstep plot.
        k = [1:r; 1:r]; k = k(:)';
        x = theta(k(2:end));
        y = prob(k(1:end-1));
        % Ensure endpoints plotted correctly.
        if x(1) >= 1 + tol, x = [1 x(1) x]; y = [0 0 y]; end
        if x(end) < th_max - tol, x = [x th_max]; y = [y y(end)]; end
        last_pos = min(find(x > th_max));
        if(isempty(last_pos)), last_pos = length(x); end
        x = x(1:last_pos);
        y = y(1:last_pos);
        x(last_pos) = th_max;
        h(j) = plot(x,y,strings{j},'Color',colors(j,:), 'Linewidth', linewidth); hold on
    end
    hold off
    xlim([1 th_max])