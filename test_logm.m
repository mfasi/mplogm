warning('off');

if (~exist('paused', 'var'))
    paused = true;
end

M = ids_max - ids_min + 1;
P = length(p_vec);
condest1 = mp(zeros(1, M));

smax = 100;
mmax_diag = 200;
mmax_tayl = 400;

condest1u = mp(zeros(P, M));
fwbound_logmct = mp(zeros(P, M));
fwbound_logagm = mp(zeros(P, M));
fwbound_logt = mp(zeros(P, M));
fwbound_logft = mp(zeros(P, M));
fwbound_logp = mp(zeros(P, M));
fwbound_logfp = mp(zeros(P, M));

for j = 1:P
    p = 2^p_vec(j);
    pp = floor(p * 2);

    ymin = ymins(j);
    ymax = ymaxs(j);
	print_legend = false;
    if (j == 1)
        print_legend = true;
    end

    mp.Digits(p);

    s_logt = zeros(1, M);
    s_logft = zeros(1, M);
    s_logp = zeros(1, M);
    s_logfp = zeros(1, M);

    m_logt = zeros(1, M);
    m_logft = zeros(1, M);
    m_logp = zeros(1, M);
    m_logfp = zeros(1, M);

    n_mats = ids_max - ids_min + 1;

    for k = ids_min : ids_max

        mp.Digits(p);
        epsilon = mp('eps');

        [A, set, id, nmats] = logm_testmats(k, n);

        % Compute reference solution.
        mp.Digits(pp);
        [exact, s, m] = logm_mp(mp(A), 'approx', 'diagonal');
        mp.Digits(p);

        norm_exact = norm(exact, 1);
        normA = norm(mp(A), 1);

         % Check backward bound of reference solution.
         mp.Digits(pp);
         bck_error_exact = norm(expm(exact) - A, 1) / normA;
         mp.Digits(p);
         log_mct = logm(mp(A));

         log_agm = logm_agm_mp(mp(A));
         try
             log_agm = logm_agm_mp(mp(A));
             fwbound_logagm(j, k) = norm(log_agm - exact, 1) / norm_exact;
         catch
             fwbound_logagm(j, k) = ymax;
         end

         [log_new_taylor, s_logt(k), m_logt(k)] = logm_mp(mp(A),...
             'approx', 'taylor',...
             'maxsqrt', smax,...
             'maxdegree', mmax_tayl);

         [log_new_full_taylor, s_logft(k), m_logft(k)] = logm_mp_full(mp(A),...
             'approx', 'taylor',...
             'maxsqrt', smax,...
             'maxdegree', mmax_tayl);

         [log_new_pade, s_logp(k), m_logp(k)] = logm_mp(mp(A),...
             'approx', 'diagonal',...
             'maxsqrt', smax,...
             'maxdegree', mmax_diag);

         [log_new_full_pade, s_logfp(k), m_logfp(k)] = logm_mp_full(mp(A),...
             'approx', 'diagonal',...
             'maxsqrt', smax,...
             'maxdegree', mmax_diag);

         if (j == 1)
             old_d = mp.Digits();
             mp.Digits(16);
             condest1(k) = funm_condest1(mp(A), @logm_mp);
             mp.Digits(p);
         end
         condest1u(j, k) = condest1(k) * epsilon;

         mp.Digits(pp);
         fwbound_logmct(j, k) = norm(log_mct - exact, 1) / norm_exact;
         fwbound_logt(j, k) = norm(log_new_taylor - exact, 1) / norm_exact;
         fwbound_logft(j, k) = norm(log_new_full_taylor - exact, 1) / norm_exact;
         fwbound_logp(j, k) = norm(log_new_pade - exact, 1) / norm_exact;
         fwbound_logfp(j, k) = norm(log_new_full_pade - exact, 1) / norm_exact;

    end
end

save('test_logm_data.mat', 'condest1u', 'fwbound_logmct', 'fwbound_logagm',...
    'fwbound_logt', 'fwbound_logft', 'fwbound_logp', 'fwbound_logfp');

%% Plot results
for j = 1:P

    print_legend = false;
    if (j == 1)
        print_legend = true;
    end

    p = 2^p_vec(j);
    mp.Digits(p);
    perfprof_norm = @(x)(max(x, x * (1 - 5e-2) + 5e-2*mp('eps')/2));

    % performance profiles
    figure(2)
    clf
    T = [perfprof_norm(fwbound_logt(j,:)'), perfprof_norm(fwbound_logp(j,:)'),...
        perfprof_norm(fwbound_logft(j,:)'), perfprof_norm(fwbound_logfp(j,:)'),...
        perfprof_norm(fwbound_logmct(j,:)'), perfprof_norm(fwbound_logagm(j,:)')];
    Tcolors = [color_logt; color_logp; color_logft; color_logfp; color_logmct; color_logagm];
    Tstyles = {ls_logt, ls_logp, ls_logft, ls_logfp, ls_logmct, ls_logagm};
    perfprof(T, 25, Tcolors, Tstyles, 1);
    axis([1 25 0 1])
    if(print_legend)
        legend('logt', 'logp', 'logft', 'logfp', 'logmct', 'logagm', 'Location', 'SE');
    end
    xlabel('theta');
    filedir = 'figs';
    filename = sprintf('%s/logm_mp_perfprof_%03d.tikz', filedir, p);
    matlab2tikz(filename, 'showInfo', false)

    ymin = ymins(j);
    ymax = ymaxs(j);

    fwbound_logmct(j, fwbound_logmct(j,:) > ymax) = ymax;
    fwbound_logagm(j, fwbound_logagm(j,:) > ymax) = ymax;
    fwbound_logt(j, fwbound_logt(j,:) > ymax) = ymax;
    fwbound_logft(j, fwbound_logft(j,:) > ymax) = ymax;
    fwbound_logp(j, fwbound_logp(j,:) > ymax) = ymax;
    fwbound_logfp(j, fwbound_logfp(j,:) > ymax) = ymax;

    fwbound_logmct(j, fwbound_logmct(j,:) < ymin) = ymin;
    fwbound_logagm(j, fwbound_logagm(j,:) < ymin) = ymin;
    fwbound_logt(j, fwbound_logt(j,:) < ymin) = ymin;
    fwbound_logft(j, fwbound_logft(j,:) < ymin) = ymin;
    fwbound_logp(j, fwbound_logp(j,:) < ymin) = ymin;
    fwbound_logfp(j, fwbound_logfp(j,:) < ymin) = ymin;

    % forward error
    figure(1)
    clf;
    [s, perm] = sort(condest1u(1,:), 'descend');
    mats = ids_min:1:ids_max;
    mp_semilogy(mats, condest1u(j, perm), [], ls_cond, 'Color', color_cond, 'MarkerSize', msize, 'Linewidth', lw);
    hold on;
    mp_semilogy(mats, fwbound_logmct(j, perm), [], 0, marker_logmct, 'Color', color_logmct, 'MarkerSize', msize);
    mp_semilogy(mats, fwbound_logagm(j, perm), [], 0, marker_logagm, 'Color', color_logagm, 'MarkerSize', msize);
    mp_semilogy(mats, fwbound_logt(j, perm), [], 0, marker_logt, 'Color', color_logt, 'MarkerSize', msize);
    mp_semilogy(mats, fwbound_logft(j, perm), [], 0, marker_logft, 'Color', color_logft, 'MarkerSize', msize);
    mp_semilogy(mats, fwbound_logp(j, perm), [], 0, marker_logp, 'Color', color_logp, 'MarkerSize', msize);
    mp_semilogy(mats, fwbound_logfp(j, perm), [0, n_mats + 1, ymin, ymax], 4, marker_logfp, 'Color', color_logfp, 'MarkerSize', msize);
    mp.Digits(p);
    if(print_legend)
        legend('klogAu', 'logmct', 'logagm', 'logt', 'logft', 'logp', 'logfp');
    end
    xlabel('phantomtheta');
    filename = sprintf('%s/logm_mp_accuracy_%03d.tikz', filedir, p);
    matlab2tikz(filename, 'showInfo', false,...
        'extraTikzpictureOptions', {'trim axis left', 'trim axis right'})

    if (paused)
        pause;
    end

end
