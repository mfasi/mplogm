warning('off');

p = 16;
pp = 2 * p;

compute_condest = true;

n_mats = ids_max - ids_min + 1;
Y = zeros(n_mats, 2, 2);
fwbound_logm_pade = zeros(1, n_mats);
fwbound_logm_tf_pade = zeros(1, n_mats);
fwbound_logm = zeros(1, n_mats);
bwbound_logm_pade = zeros(1, n_mats);
bwbound_logm_tf_pade = zeros(1, n_mats);
bwbound_logm = zeros(1, n_mats);
condest1u = zeros(1, n_mats);

group_ids = cell(n_mats, 1);

for i = ids_min : ids_max

    [A, set, id, nmats] = logm_testmats(i, n);

    mp.Digits(p);
    [log_logm_pade, s_logm_pade, m_logm_pade] = logm_mp(double(A),...
        'approx', 'diagonal',...
        'maxsqrt', 50,...
        'maxdegree', 100);
    [log_logm_pf_pade, s_logm_tf_pade, m_logm_tf_pade] = logm_mp_full(double(A),...
        'approx', 'diagonal',...
        'maxsqrt', 50,...
        'maxdegree', 100);
    [~, s_in_logm, s_logm, m_logm] = logm_2011(double(A));
    [log_logm] = logm(double(A));

    mp.Digits(2 * pp);
    exact = logm_mp(mp(A),...
        'approx', 'taylor');

    fwbound_logm_pade(i) = double(norm(log_logm_pade - exact, 1) / norm(exact, 1));
    fwbound_logm_tf_pade(i) = double(norm(log_logm_pf_pade - exact, 1) / norm(exact, 1));
    fwbound_logm(i) = double(norm(log_logm - exact, 1) / norm(exact, 1));

    if (compute_condest)
        mp.Digits(p);
        condest1u(i) = funm_condest1(mp(A), @logm_mp, @logm_mp) * eps();
    end

    Y(i, 1, 2) = s_logm_pade;
    Y(i, 1, 3) = m_logm_pade;
    Y(i, 2, 2) = s_logm;
    Y(i, 2, 3) = m_logm;

    group_ids{i} = i;

end

% performance profile
perfprof_norm = @(x)(max(x, x * (1 - 5e-2) + 5e-2*eps/2));
figure(10)
clf
T = [perfprof_norm(fwbound_logm)',...
     perfprof_norm(fwbound_logm_pade)',...
     perfprof_norm(fwbound_logm_tf_pade)'];
Tcolors = [color_logm; color_logp; color_logfp];
Tstyles = {ls_logm; ls_logp; ls_logfp};
perfprof(T, 25, Tcolors, Tstyles, 1);
axis([1 25 0 1])
legend('logm', 'logp', 'logfp', 'Location', 'SE');
xlabel('theta');

% save .png
saveas(gcf, 'pngfigs/bck_fwd_perfprof', 'png');

% save .tikz
matlab2tikz('figs/bck_fwd_perfprof.tikz', 'showInfo', false,...
        'extraTikzpictureOptions', {'trim axis left', 'trim axis right'});

min_y = 1e-18;
[~, perm] = sort(condest1u, 'descend');
fwbound_logm_pade(fwbound_logm_pade < min_y) = min_y;
fwbound_logm_tf_pade(fwbound_logm_tf_pade < min_y) = min_y;
fwbound_logm(fwbound_logm < min_y) = min_y;

% conditioning and forward error
figure(1)
clf
if (compute_condest)
    p1 = semilogy(condest1u(perm), ls_cond, 'Color', color_cond, 'LineWidth', lw);
    hold on;
else
    p1 = [];
end
p2 = semilogy(fwbound_logm(perm), marker_logm, 'Color', color_logm);
hold on;
p3 = semilogy(fwbound_logm_pade(perm), marker_logp, 'Color', color_logp);
p4 = semilogy(fwbound_logm_tf_pade(perm), marker_logfp, 'Color', color_logfp);
axis([0 ids_max+1 min_y 1])
h1 = gca;
legend([p1, p2, p3, p4], 'klogAu',  'logm', 'logp', 'logfp');

% save .png
if png_out
    saveas(gcf, 'pngfigs/bck_fwd_plot', 'png');
end

% save .tikz
if tikz_out
    matlab2tikz('figs/bck_fwd_plot.tikz', 'showInfo', false,...
        'extraTikzpictureOptions', {'trim axis left', 'trim axis right'});
end

figure(3)
clf
tot_vec = (Y(perm, 1, 2) + Y(perm, 1, 3))./(Y(perm, 2, 2) + Y(perm, 2, 3));
tot_vec(isnan(tot_vec)) = 1; % from 0 / 0
sqrt_vec = Y(perm, 1, 2)./Y(perm, 2, 2);
sqrt_vec(isnan(sqrt_vec)) = 1; % from 0 / 0
plot(tot_vec, 'bx', 'MarkerSize', msize);
hold on
axis([0 ids_max+1 round((min(tot_vec) - 0.05) * 20)/20 round((max(tot_vec) + 0.025) * 20)/20]);
xlabel('phantomtheta');

% save .png
if png_out
    saveas(gcf, 'pngfigs/bck_fwd_cost', 'png');
end

% save .tikz
if tikz_out
    matlab2tikz('figs/bck_fwd_cost.tikz', 'showInfo', false,...
        'extraTikzpictureOptions', {'trim axis left', 'trim axis right'});
end

% Histogram (not published)
figure(4)
clf
p4 = plotBarStackGroups_bck_fwd(Y(perm,:,:), group_ids(perm));
axis([0 ids_max+1 0 60]);

% Estimate overscaling on data set
if false
    sqrts = Y(:,:,2);
    sqrts(sqrts == 0) = 10000;
    sqrts(:,2) = min(sqrts(:,2), sqrts(:,1));
    ss = (sqrts(:,1) - sqrts(:,2)) ./ sqrts(:,1);
    nnz_ss = ss(find(ss));
    fprintf('Overscaling on %.2f%% of the data set.\n',...
        nnz(ss) / length(sqrts(:,1)) * 100);

    fprintf('Up to %.2f%% less square roots.\n',...
        max(nnz_ss) * 100);
    fprintf('On average %.2f%% less square roots.\n',...
        mean(nnz_ss) * 100);
    fprintf('At least %.2f%% less square roots.\n',...
        min(nnz_ss) * 100);
end