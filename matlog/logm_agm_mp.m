function X = logm_agm_mp(A)

n = size(A, 1);
epsilon = myeps(class(A));

vareps = sqrt(epsilon) / norm(A,'fro');

X = -mp('pi') / 2 * inv(agm_lt(vareps * A, epsilon));

X(1:n+1:n^2) = X(1:n+1:n^2) + log(4 / vareps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             SUBFUNCTIONS                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function X = agm_lt(A, epsilon)
        % This implements the Legendre variant for the computation of the arithmetic
        % geometric mean of a matrix A such that ||A||^2 <= u(unit roundoff). It will
        % provide innacurate results if the previous condition is not met because it
        % approximates the inverse (I+A)^-1 with I-A

        n = size(A, 1);
        I = eye(n);
        k = 1;
        [S, iterDB(k)] = sqrtm_dbp_mod(A,epsilon);
        P = 2 * S * (I - A);
        Q = (I + A) / 4 * (I + P);
        curr_err = norm(I - P, 'fro');
        maxiter = 40;

        while (k <= maxiter && curr_err > 4 * myeps(class(A)))
            k = k + 1;
            [S, iterDB(k)] = sqrtm_dbp_mod(P, epsilon);
            P = 2 * S / (I + P);
            Q = Q * (I + P) / 2;
            curr_err = norm(I-P,'fro');
        end

        X = Q;

    end

    function [X,k] = sqrtm_dbp_mod(A, epsilon)
        % This is a modification of sqrtm_dbp from MFToolbox that avoids overflow
        % and underflow in the computation of the scaling factor
        % k is the number of iterations

        n = size(A,1);
        I = eye(n);

        tol = mft_tolerance(A, epsilon);
        X = A;
        M = A;
        maxit = 25;
        scale = 1;

        for k = 1:maxit

            if scale == 1
                [~, U] = lu(M);
                g = prod(abs(diag(U)).^(-1/(2*n)));
                X = g * X;
                M = g^2 * M;
            end

            Xold = X;
            invM = inv(M);

            X = X * (I + invM)/2;
            M = 0.5 * (I + (M + invM)/2);

            Mres = norm(M - I,'fro');

            reldiff = norm(X - Xold,'fro') / norm(X,'fro');
            if reldiff < 1e-2
                scale = 0;
            end

            if Mres <= tol
                return;
            end

        end
        error('Not converged after %2.0f iterations', maxit)
    end

    function e = myeps(curr_class)
        %COMP_EPS    Machine epsilon.
        %   COMP_EPS(CLASS) computes the machine epsilon of class CLASS.
        if(strcmp(curr_class, 'mp'))
            e = mp('eps');
        else
            e = eps(curr_class);
        end
    end

    function tol = mft_tolerance(A, epsilon)
        %MFT_TOLERANCE   Convergence tolerance for matrix iterations.
        %   TOL = MFT_TOLERANCE(A) returns a convergence tolerance to use in
        %   the matrix iterations in the Matrix Function Toolbox applied to the
        %   matrix A.  All functions in the toolbox call this function to set
        %   the convergence tolerance.

        tol = sqrt(length(A))*epsilon/2;
    end

end