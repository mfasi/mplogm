sizes = [10, 20, 50, 100, 200, 500];
nsizes = length(sizes);

p = 512;
mp.Digits(p);

filename = sprintf('tabs/table_profile_td_%04d.tex', p);
fileid_td = fopen(filename,'w');
fprintf(fileid_td, ['\\begin{tabularx}{\\textwidth}',...
    '{@{\\extracolsep{\\fill}}rr|rrrrrrr|rrrrrrr}\n']);
fprintf(fileid_td, '\\toprule\n');
fprintf(fileid_td, ['\\multicolumn{2}{c|}{} & ',...
    '\\multicolumn{7}{c|}{\\logt} & ',...
    '\\multicolumn{7}{c}{\\logp} \\\\\n']);
fprintf(fileid_td, [' & $n$ & ',...
    '$s$ & $m$ &  $T_{sch}$ & $T_{sqrt}$ & $T_{bnd}$ & $T_{eval}$ & $T_{tot}$ & ',...
    '$s$ & $m$ &  $T_{sch}$ & $T_{sqrt}$ & $T_{bnd}$ & $T_{eval}$ & $T_{tot}$ \\\\\n']);
fprintf(fileid_td, '\\midrule\n');

mats_id = {'\verb|A|', '\verb|B|', '\verb|C|'};

nmat = 3;
for j = nmat %1:nmat

    for i = nsizes %1:nsizes

        n = sizes(i);

        switch j
            case 1
                A = expm(gallery('randsvd', n));
            case 2
                A = expm(gallery('chebvand', n));
            case 3
                A = expm(gallery('chow', n));
        end

        [X, s_t, m_t, time_t] = logm_mp(mp(A),...
            'approx', 'taylor',...
            'epsilon', mp('eps'),...
            'maxsqrt', 80,...
            'maxdegree', 40,...
            'timing', true);

        [X, s_d, m_d, time_d] = logm_mp(mp(A),...
            'approx', 'diagonal',...
            'epsilon', mp('eps'),...
            'maxsqrt', 80,...
            'maxdegree', 40,...
            'timing', true);

        s1 = sum(time_t)/100;
        s2 = sum(time_d)/100;

        if (i == 1)
            id_field = mats_id{j};
        else
            id_field = '            ';
        end
        fprintf(fileid_td,...
            '%s & %3d & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f \\\\\n',...
            id_field, sizes(i),...
            sum(s_t), m_t, time_t(1)/s1, sum(time_t(2:3))/s1, time_t(4)/s1 , time_t(5)/s1,s1*100,...
            sum(s_d), m_d, time_d(1)/s2, sum(time_d(2:3))/s2, time_d(4)/s2, time_d(5)/s2,s2*100);

        fprintf('%s & %3d & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f & %2d & %2d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f \\\\\n',...
            id_field, sizes(i),...
            sum(s_t), m_t, time_t(1)/s1, sum(time_t(2:3))/s1, time_t(4)/s1 , time_t(5)/s1,s1*100,...
            sum(s_d), m_d, time_d(1)/s2, sum(time_d(2:3))/s2, time_d(4)/s2, time_d(5)/s2,s2*100);
    end
    if (i == nsizes && j ~= nmat)
        fprintf(fileid_td, '\\midrule\n');
    end
end

fprintf(fileid_td, '\\bottomrule\n');
fprintf(fileid_td, '\\end{tabularx}');