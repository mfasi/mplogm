for id = [1, 2, 3]

    n = 10;
    switch id
        case 1
            A = expm([4 2 0; 1 4 1; 1 1 4]);
        case 2
            A = expm(gallery('dramadah', n));
        case 3
            A = expm(gallery('toeppen', n));
    end

    precisions = 2.^(1:1:12); % could try 1:1:11
    n_prec = length(precisions);

    smax = 80;
    mmax = 200;

    fwbound_logpr = mp(zeros(1, n_prec));
    s_logpr = mp(zeros(1, n_prec));
    m_logpr = mp(zeros(1, n_prec));
    fwbound_logpa = mp(zeros(1, n_prec));
    s_logpa = mp(zeros(1, n_prec));
    m_logpa = mp(zeros(1, n_prec));

    fwbound_logtr = mp(zeros(1, n_prec));
    s_logtr = mp(zeros(1, n_prec));
    m_logtr = mp(zeros(1, n_prec));
    fwbound_logta = mp(zeros(1, n_prec));
    s_logta = mp(zeros(1, n_prec));
    m_logta = mp(zeros(1, n_prec));

    for i = 1:n_prec

        p = precisions(i);

        % Compute reference solution
        mp.Digits(2 * p);
        [exact, s, m] = logm_mp(mp(A),...
            'approx', 'diagonal',...
            'epsilon', mp('eps'),...
            'maxsqrt', smax,...
            'maxdegree', mmax);
        mp.Digits(p);
        norm_exact = norm(exact, 1);

        [logp_rel, s_logpr(i), m_logpr(i)] = logm_mp(mp(A),...
            'approx', 'diagonal',...
            'epsilon', mp('eps'),...
            'maxsqrt', smax,...
            'maxdegree', mmax);

        [logp_abs, s_logpa(i), m_logpa(i)] = logm_mp_abs(mp(A),...
            'approx', 'diag',...
            'epsilon', mp('eps'),...
            'maxsqrt', smax,...
            'maxdegree', mmax);

        [logt_rel, s_logtr(i), m_logtr(i)] = logm_mp(mp(A),...
            'approx', 'taylor',...
            'epsilon', mp('eps'),...
            'maxsqrt', smax,...
            'maxdegree', mmax);

        [logt_abs, s_logta(i), m_logta(i)] = logm_mp_abs(mp(A),...
            'approx', 'taylor',...
            'epsilon', mp('eps'),...
            'maxsqrt', smax,...
            'maxdegree', mmax);
        fwbound_logpr(i) = (norm(logp_rel - exact, 1) / norm_exact) / mp('eps');
        fwbound_logpa(i) = (norm(logp_abs - exact, 1) / norm_exact) / mp('eps');
        fwbound_logtr(i) = (norm(logt_rel - exact, 1) / norm_exact) / mp('eps');
        fwbound_logta(i) = (norm(logt_abs - exact, 1) / norm_exact) / mp('eps');
    end

    figure(1)
    clf

    loglog(precisions, fwbound_logpa, 'r-', 'Linewidth', lw);
    hold on
    loglog(precisions, fwbound_logta, 'm--', 'Linewidth', lw);
    loglog(precisions, fwbound_logpr, 'g-.', 'Linewidth', lw);
    loglog(precisions, fwbound_logtr, 'b:', 'Linewidth', lw);
    loc_const = double(max([fwbound_logpa, fwbound_logpr]));
    xlabel('d');
    axis([min(precisions), max(precisions), 1e-1, 1e15]);

    legend('logpa', 'logta', 'logpr', 'logtr', 'Location', 'NW');

    % save .png
    filename = sprintf('pngfigs/abs_rel_def_%02d', id);
    saveas(gcf, filename, 'png');

    % save .tikz
    filedir = 'figs';
    filename = sprintf('%s/abs_rel_%02d.tikz', filedir, id);
    matlab2tikz(filename, 'showInfo', false,...
        'extraTikzpictureOptions', {'trim axis left', 'trim axis right'});

end
