function optimize_error_ratio()
    n = 5;
    A0 = randn(n) + n * eye(n);
    f_err_fwd_bck = @(fbound_fwd, sm_fwd, fbound_bck, sm_bck)...
        (fbound_fwd/fbound_bck);
    f_err_bck_fwd = @(fbound_fwd, sm_fwd, fbound_bck, sm_bck)...
        (fbound_bck/fbound_fwd);

    [x, f_err_fwd_bck_max, nf] =...
        mdsmax(@(X)(optimize_subroutine(X, f_err_fwd_bck)), A0);
    [x, f_err_bck_fwd_max, nf] =...
        mdsmax(@(X)(optimize_subroutine(X, f_err_bck_fwd)), A0);

    f_cost_fwd_bck = @(fbound_fwd, sm_fwd, fbound_bck, sm_bck)...
        (sm_fwd/sm_fwd);
    f_cost_bck_fwd = @(fbound_fwd, sm_fwd, fbound_bck, sm_bck)...
        (sm_fwd/sm_fwd);

    [x, f_cost_fwd_bck_max, nf] =...
        mdsmax(@(X)(optimize_subroutine(X, f_cost_fwd_bck)), A0);
    [x, f_cost_bck_fwd_max, nf] =...
        mdsmax(@(X)(optimize_subroutine(X, f_cost_bck_fwd)), A0);

    f_err_fwd_bck_max
    f_err_bck_fwd_max
    f_cost_fwd_bck_max
    f_cost_bck_fwd_max


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function er = optimize_subroutine(A, f)
        % A: The matrix
        % f: A function that takes A

        eigs_v = eig(A);
        if (any(imag(eigs_v) == 0 & eigs_v <= 0))
            er = -inf;
        else
            mp.Digits(16);
            [log_fwd, s_fwd, m_fwd] = logm_mp(double(A),...
                'approx', 'diagonal',...
                'epsilon', eps(),...
                'maxsqrt', 50,...
                'maxdegree', 100);
            [log_bck, s_bck, m_bck] = logm_2011(double(A), 50);

            mp.Digits(32);

            exact = logm(mp(A));
            fwbound_fwd = norm(log_fwd - exact) / norm(exact);
            fwbound_bck = norm(log_bck - exact) / norm(exact);

            er = f(fwbound_fwd,s_fwd + m_fwd,...
                fwbound_bck, s_bck + m_bck);
        end
        end
end
