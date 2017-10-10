function mp_semilogy(x, y, axis_lim, n_yticks, varargin)
% MP_SEMILOGY Handle values smaller than 1e-300 with SEMILOGY.

% Copyright (C) 2016 Massimiliano Fasi

    plot(x, log10(y), varargin{:})
    if (isempty(axis_lim))
        axis([min(x)-1, max(x)+1, double(min(log10(y(y>0))))-1, double(max(log10(y(y>0))))])
    else
        axis_lim(3:4) = log10(axis_lim(3:4));
        axis_lim = double(axis_lim);
        axis(axis_lim);
        if (n_yticks > 0)
            step = (axis_lim(4) - axis_lim(3)) / (n_yticks - 1);
            set(gca,'Ytick', axis_lim(3):step:axis_lim(4))
        end
    end
    ax = gca;
    y = get(ax, 'YTick');
    strings = cell(length(y), 1);
    for i = 1 : length(y)
        if(ceil(y(i)) == y(i))
            strings{i} = sprintf('10^{%d}', y(i));
        else
            strings{i} = sprintf('10^{%.1f}', y(i));
        end
    end
    set(ax,'YTickLabel',strings)
end
