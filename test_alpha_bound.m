mp.Digits(512);

A1 = mp([0.01 0.95; 0.00 0.04]);

m_vec = 1:1:20;
m_len = length(m_vec);

alpha_bound_diag_minus_A1 = zeros(m_len, 1);
norm_bound_diag_minus_A1 = zeros(m_len, 1);
fw_norm_taylor_A1 = zeros(m_len, 1);

alpha_bound_taylor_minus_A1 = zeros(m_len, 1);
norm_bound_taylor_minus_A1 = zeros(m_len, 1);
fw_norm_diag_A1 = zeros(m_len, 1);

for i = 1 : m_len

    m = m_vec(i);

    % taylor
    d_max = floor((1 + sqrt(4 * m + 5)) / 2);

    alpha_taylor_A1 = Inf;
    for d = 1 : d_max
        dp1 = normAm(A1, d)^(1/d);
        dp2 = normAm(A1, d+1)^(1/d+1);
        new_alpha = max(dp1, dp2);
        alpha_taylor_A1 = min(alpha_taylor_A1, new_alpha);
    end
    normA1 = norm(A1, 1);

    v = mp(m:-1:1);
    p = ((-1).^(v+1))./v;
    p = [p, 0];

    fw_norm_taylor_A1(i) = norm(logm_taylor(A1, i) - logm(eye(2) + A1), 1);
    alpha_bound_taylor_minus_A1(i) =...
        abs(polyval(p, -alpha_taylor_A1) - log(1 - alpha_taylor_A1));
    norm_bound_taylor_minus_A1(i) =...
        abs(polyval(p, -normA1) - log(1 - normA1));

    % diagonal
    d_max = floor((1 + sqrt(8 * m + 5)) / 2);

    alpha_diag_A1 = Inf;
    for d = 1 : d_max
        dp1 = normAm(A1, d)^(1/d);
        dp2 = normAm(A1, d+1)^(1/d+1);
        new_alpha = max(dp1, dp2);
        alpha_diag_A1 = min(alpha_diag_A1, new_alpha);
    end

    [nodes, wts] = mp.GaussLegendre(m, 0, 1);
    S_alpha_minus_A1 = 0;
    S_norm_minus_A1 = 0;
    S_fw_norm_A1 = zeros(2, 2);
    for j = 1:m
        S_alpha_minus_A1 = S_alpha_minus_A1 + wts(j) * (-alpha_diag_A1 / (1 - nodes(j) * alpha_diag_A1));
        S_norm_minus_A1 = S_norm_minus_A1 + wts(j) * (-normA1 / (1 - nodes(j) * normA1));
        S_fw_norm_A1 = S_fw_norm_A1 + wts(j) * (A1 / (eye(2) + nodes(j) * A1));
    end

    fw_norm_diag_A1(i) = norm(S_fw_norm_A1 - logm(eye(2) + A1), 1);
    alpha_bound_diag_minus_A1(i) =...
        abs(S_alpha_minus_A1 - log(1 - alpha_diag_A1));
    norm_bound_diag_minus_A1(i) =...
        abs(S_norm_minus_A1 - log(1 - normA1));

end

x = [m_vec(1), m_vec(end)];

color_normminus = [0 0 0.5];
color_alphaminus = [0 0.7 0];
color_fw = [0.8 0 0];

ls_normminus = ':';
ls_alphaminus = '-.';
ls_fw = '-';

% taylor
figure(1)
clf
semilogy(m_vec, norm_bound_taylor_minus_A1,...
    ls_normminus, 'Color', color_normminus, 'Linewidth', lw);
hold on;
semilogy(m_vec, alpha_bound_taylor_minus_A1,...
    ls_alphaminus, 'Color', color_alphaminus, 'Linewidth', lw);
semilogy(m_vec, fw_norm_taylor_A1,...
    ls_fw, 'Color', color_fw, 'Linewidth', lw);
axis([1 m_vec(end) 1e-98 1e2])
legend('normminust', 'alphaminust', 'forwarderrort',...
    'Location', 'SouthWest')
xlabel('m');

% save .png
if png_out
    saveas(gcf, 'pngfigs/alpha_norm_taylor', 'png');
end

% save .tikz
if tikz_out
    matlab2tikz('figs/alpha_norm_taylor.tikz', 'showInfo', false,...
        'extraTikzpictureOptions', {'trim axis left', 'trim axis right'});
end

% diagonal
figure(2)
clf
semilogy(m_vec, norm_bound_diag_minus_A1,...
    ls_normminus, 'Color', color_normminus, 'Linewidth', lw);
hold on;
semilogy(m_vec, alpha_bound_diag_minus_A1,...
    ls_alphaminus, 'Color', color_alphaminus, 'Linewidth', lw);
semilogy(m_vec, fw_norm_diag_A1,...
    ls_fw, 'Color', color_fw, 'Linewidth', lw);
axis([1 m_vec(end) 1e-98 1e2])
legend('normminusp', 'alphaminusp', 'forwarderrorp',...
    'Location', 'SouthWest')
xlabel('m');

% save .png
if png_out
    saveas(gcf, 'pngfigs/alpha_norm_diagonal', 'png');
end

% save .tikz
if tikz_out
    matlab2tikz('figs/alpha_norm_diagonal.tikz', 'showInfo', false,...
        'extraTikzpictureOptions', {'trim axis left', 'trim axis right'});
end

function S = logm_taylor(A,m)
%LOGM_TAYLOR   Taylor approximation to matrix logarithm.
%   LOGM_TAYLOR(A,M) computes the Taylor approximant of degree M to
%   LOG(EYE(SIZE(A))+A) using Paterson-Stockmeyer algorithm.

    v_0m = cast(0:1:m, class(A));
    c = (1./v_0m').*(-1).^(v_0m' + 1);
    c(1) = zeros(1, class(A));

    ss = ceil(sqrt(m));
    r = floor(m/ss);

    Xpow = cell(ss+1, 1);
    id_n_loc = eye(size(A,1), class(A));
    Xpow{1} = id_n_loc;
    Xpow{2} = A;
    for i=3:ss+1
        Xpow{i} = A*Xpow{i-1};
    end

    B = cell(r+1);
    for kk = 0 : r-1
        tt = c(ss*kk+1) * id_n_loc;
        for j=1:ss-1
            tt = tt + c(ss*kk+j+1)*Xpow{j+1};
        end
        B{kk+1} = tt;
    end

    B{r+1} = c(m+1)*Xpow{m-ss*r+1};
    for j=m-1:-1:ss*r
        if j == ss*r
           B{r+1} = B{r+1} + c(ss*r+1)*id_n_loc;
        else
           B{r+1} = B{r+1} + c(j+1)*Xpow{m-ss*r-(m-j)+1};
        end
    end

    As = Xpow{ss+1};
    S = zeros(size(A,1),class(A));
    for kk=r:-1:0
        S = S * As + B{kk+1};
    end
end
