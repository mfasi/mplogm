function [varargout] = logm_mp(A, varargin)
%LOGM_MP  Scur algorithm for multiprecision matrix logarithm and derivative.
%   This code requires the Advanpix Multiprecision Computing Toolbox
%   (see www.advanpix.com).
%
%   [X, S, M] = logm_mp(A) computes the principal matrix logarithm X of the
%   matrix A.
%
%   [Y, S, M] = logm_mp(A,E) computes the Frechet derivative Y of the
%   matrix logarithm of A in the direction E.
%
%   [X, Y, S, M] = logm_mp(A,E) computes both the principal matrix
%   logarithm X of the matrix A and the Frechet derivative Y of A in the
%   direction E.
%
%   The output parameters S and M are the number of square roots taken by
%   the algorithm and the degree of the Pade approximation, respectively.
%
%   [...] = logm_mp(...,'precision',DIGITS) specifies the number of digits
%   to be used in the computation. Default is mp.Digits() if A
%   is of class mp. The computation is performed in single or double
%   arithmetic if A is of class single or double, respectively, and
%   the precision is not specified.
%
%   [...] = logm_mp(...,'epsilon',EPSILON) specifies the tolerance to be
%   used to evaluate the Pade approximations. Default is machine epsilon
%   of the precision of A if A is of class 'single' or 'double', mp.eps()
%   if A is of class 'mp'.
%
%   [...] = logm_mp(...,'maxsqrtm',MAXSQRTM) specifies the maximum number
%   of square roots to be taken. Default is 100.
%
%   [...] = logm_mp(...,'maxdegree',MAXDEGREE) specifies the maximum
%   degree of the Pade approximants. Default is 100.
%
%   [...] = logm_mp_full(...,'approx', KIND), where KIND='diagonal' or
%   KIND='taylor', specifies the kind of Pade approximant to use. Default
%   is diagonal. The efficient evaluation of truncated Taylor series
%   requires storing several additional matrices having the same size as A.
%
%   Reference: Algorithm 4.1 in
%   M. Fasi and N. J. Higham, Multiprecision Algorithms for Computing the
%   Matrix Logarithm. Technical Report 2017.16, Manchester Institute for
%   Mathematical Sciences, The University of Manchester, UK, May 2017.

% Copyright (c) 2016-2017, Massimiliano Fasi
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
% OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%% Check whether the Advanpix MCT toolbox is available.
if ~check_tb_available('Advanpix Multiprecision Computing Toolbox')
  err_msg = ['The function LOGM_MP requires the Advanpix '...
    'Multiprecision Computing Toolbox,\nwhich does not appear '...
    'to be installed on this system.\nSee www.advanpix.com.'];
  error('logm_mp:missing_toolbox', err_msg);
end

%% Parse and validate input.
p = inputParser;
addRequired(p, 'A', @(x)(ismatrix(x) && isnumeric(x)));
addOptional(p, 'E', [], @(x)(ismatrix(x) && isnumeric(x)));
addParameter(p, 'precision', [], @(x)(x == round(x) && x > 0))
addParameter(p, 'epsilon', [], @(x)(x > 0 && x < 1));
addParameter(p, 'maxsqrt', 100, @(x)(x == round(x) && x > 0));
addParameter(p, 'maxdegree', 100, @(x)(x == round(x) && x > 0));
addParameter(p, 'approx', 'diagonal',...
  @(x)(ischar(x) && strcmp(x, 'taylor') || strcmp(x, 'diagonal')));
addParameter(p, 'timing', false, @(x)(x == true || x == false));

parse(p,A,varargin{:});

A = p.Results.A;
E = p.Results.E;
digits = p.Results.precision;
epsilon = p.Results.epsilon;
maxsqrt = p.Results.maxsqrt;
maxdegree = p.Results.maxdegree;
timing = p.Results.timing;

if strcmp(p.Results.approx, 'taylor')
  usetaylor = true;
  eval_pade_error = @scalar_error_taylor;
  alpha_pade = @(A, m)(alpha(A,m,0));
else
  usetaylor = false;
  eval_pade_error = @scalar_error_diagonal;
  alpha_pade = @(A, m)(alpha(A,m,m));
end

min_pade_deg = @(A, epsilon, maxd)...
  (find_min_deg(A, epsilon, maxd, eval_pade_error, alpha_pade));

[n1, n] = size(A);
if n ~= n1
  error('The matrix ``A'' must be square.');
end

%% Parse and validate output.
if nargout > 4
  error('This function returns at most four values.');
end
if nargout < 4 && timing
  error('With (''timing'', true), the function returns four values.');
end
compute_frechet = false;
if ~isempty(E)
  compute_frechet = true;
  [mE,nE] = size(E);
  if nE ~= mE
    error('The matrix ``E'' must be square.');
  end
  if nE ~= n
    error('The matrices ``A'' and ``E'' must have the same size.');
  end
end
compute_logm = ~(compute_frechet &&...
  (nargout == 3 || nargout == 1 || nargout == 0));

%% Determine whether the computation will be single, double or mp.
current_digits = mp.Digits();
if ~isempty(digits) % working precision specified
  mp.Digits(digits)
  A = mp(A);
else
  if isa(A, 'double')
    digits = 16;
  elseif isa(A, 'single')
    digits = 8;
  else
    digits = mp.Digits();
  end
  mp.Digits(digits);
end
if isempty(epsilon)
  epsilon = myeps(class(A));
end

%% Form complex Schur form (if A not already upper triangular).
id_n = eye(n, class(A));
schur_time = tic;
schurInput = istriu(A);
if schurInput
  T = A;
  Q = id_n;
else
  [Q, T] = schur(A, 'complex');
  if compute_frechet
    E = Q' * E * Q;
  end
end
T0 = T;
time(1) = toc(schur_time);

use_ccc = false;
if any(imag(diag(T)) == 0 & real(diag(T)) <= 0)
  use_ccc = true;
  warning(['Principal logarithm not defined for matrices with\n',...
    'nonpositive real eigenvalues, the solution conforms\n',...
    'to the rule of counter-clockwise continuity.']);
end

s = 0;
if isdiag(T) && ~compute_frechet      % Check if T is diagonal.
    d = diag(T);
    logd = log(d);
    X = (Q.*logd.')*Q';
    if isreal(logd)
        X = (X+X')/2;
    end
    m = 0;
else
    %% Compute initial square roots.
    % Take square root until a Pade approximant of degree maxdegree is
    % enough to achieve the required tolerance.
    init_sqrtm = tic;
    Tm1 = T - id_n;
    eps_normTm1 = epsilon * norm(Tm1, 1);
    a = -alpha_pade(Tm1, maxdegree);
    while ((max(abs(diag(Tm1))) >= 1 || abs(a) > 1 ||...
            eval_pade_error(a, maxdegree) >= eps_normTm1) &&...
            s < maxsqrt)
        [T, E, s] = compute_sqrtm(T, E, s, compute_frechet);
        Tm1 = T - id_n;
        eps_normTm1 = epsilon * norm(Tm1, 1);
        a = -alpha_pade(Tm1, maxdegree);
    end
    time(2) = toc(init_sqrtm);
    
    %% Estimate good value for m.
    time(3) = 0;
    min_prob = tic;
    eps_normTm1 = epsilon * norm(Tm1, 1);
    m = min_pade_deg(T, eps_normTm1, maxdegree);
    
    old_m = m;
    count = 0;
    try_sqrtm = true;
    while (try_sqrtm && s < maxsqrt)
        if usetaylor
            newm = max(ceil((sqrt(m) - 0.5)^2 - 1), 1);
        else
            newm = max(m - 2, 1);
        end
        a = -alpha_pade(Tm1, m) / 2;
        boundabs2 = eval_pade_error(a, newm);
        
        % Check if taking sqrtm and a lower degree approximant will work.
        if boundabs2 > eps_normTm1              % no
            try_sqrtm = false;
        else                                    % yes
            sub_sqrtm = tic;
            [T, E, s] = compute_sqrtm(T, E, s, compute_frechet);
            time(3) = time(3) + toc(sub_sqrtm);
            Tm1 = T - id_n;
            eps_normTm1 = epsilon * norm(Tm1, 1);
            m = min_pade_deg(T, eps_normTm1, m);
            if m >= old_m
                count = count + 1;
            else
                count = 0;
            end
            old_m = m;
        end
    end
    time(4) = toc(min_prob) - time(3);
    
    final_comp = tic;
    %% Compute accurate superdiagonal of T^(1/2^s).
    for l = 1:n-1
        T(l,l+1) = powerm2by2(T0(l:l+1,l:l+1), 1/2^s);
    end
    
    %% Compute accurate diagonal of T^(1/2^s) - I.
    T(1:n+1:end) = sqrt_power_1(diag(T0), s);
    
    if ~usetaylor || compute_frechet
        [nodes, wts] = mp.GaussLegendre(m, 0, 1);
        nodes = cast(nodes, class(T));
        wts = cast(wts, class(T));
    end
    
    %% Compute matrix logarithm.
    if compute_logm
        if usetaylor
            X = logm_taylor(T, m);
        else
            X = zeros(size(T,1), class(T));
            for k = 1:m
                X = X + wts(k) * ((id_n + nodes(k)*T) \ T);
            end
        end
        
        X = 2^s * X;
        
        % Compute accurate diagonal and superdiagonal of log(T).
        for l = 1:n-1
            X(l:l+1, l:l+1) = logm2by2(T0(l:l+1, l:l+1));
        end
        
        X = Q * X * Q';
        if isreal(A) && ~use_ccc
            X = real(X);
        end
    end
    time(5) = toc(final_comp);
    
    %% Compute Frechet derivative.
    if compute_frechet
        L = zeros(n,n,class(T)); % The derivative
        for k = 1:m
            temp = id_n + nodes(k) * T;
            L = L + wts(k) * (temp \ E) / temp;
        end
        
        L = 2^s * Q * L * Q';
        if isreal(A)
            L = real(L);
        end
    end
end

%% Prepare output.
k = 1;
if compute_logm
  [varargout{k}] = X;
  k = k + 1;
end
if compute_frechet
  [varargout{k}] = L;
  k = k + 1;
end
if k < nargout
  [varargout{k}] = s;
  k = k + 1;
  [varargout{k}] = m;
end
if k < nargout
  k = k + 1;
  [varargout{k}] = time;
end

%% Restore mp working precision.
mp.Digits(current_digits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             SUBFUNCTIONS                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Square root.

  function [R, E, s] = compute_sqrtm(T, E, s, compute_frechet)
    if ismp(T)
      R = sqrtm_tri(T); % this is the function in mp.m
    else
      R = sqrtm(T);     % call matfun/private/sqrtm_tri.m
    end
    if compute_frechet
      if ismp(R)
        E = sylvester_tri(R,R,E);
      else
        E = sylvester(R,R,E); % call matfun/private/sylvester_tri.m
      end
    end
    s = s + 1;
  end

%% Scalar bound evaluation and Pade approximation.

  function e = scalar_error_taylor(x, m)
    %LOGM_TAYLOR   Taylor approximation to matrix logarithm.
    %   LOGM_TAYLOR(A,M) computes the Taylor approximant of degree M to
    %   LOG(EYE(SIZE(A))+A) using Paterson-Stockmeyer algorithm.
    x = mp(x);
    y = sum(-cumprod(repmat(-x, 1, m))./(1:m));
    e = abs(y - log(1 + x));
  end

  function e = scalar_error_diagonal(x, m)
    %LOG_DIAGONAL   Pade approximation to scalar logarithm.
    %   LOGM_PADE(A,M) computes the [M/M] Pade approximant to
    %   LOG(1+A) using bottom-up evaluation of the continued
    %   fraction expansion of the Pade approximants.
    x = mp(x);
    y =  m / (4 * m - 2) * x;
    for j = m-1 : -1 : 1
      y =  (j * x) / ((j * 4 + 2) * (1 + y));
      y =  (j * x) / ((j * 4 - 2) * (1 + y));
    end
    y = x / (1 + y);
    e = abs(y - log(1 + x));
  end

  function m = find_min_deg(A, epsilon, maxd, scalar_error_pade, alpha_pade)
    %MIN_DIAGONAL_DEG   Degree of Pade approximation.
    %    MIN_DIAGONAL_DEG(T,EPSILON,MAXD) finds the minimum m such that the
    %    [m/m] Pade approximant is accurate enought to compute the matrix
    %    logarithm of T + EYE(n) with truncation error smaller than
    %    EPSILON.
    m_left = 1;
    m_right = maxd;
    found = false;
    AmI = A - eye(size(A,1), class(A));
    while (~found) % bisection
      m_middle = ceil((m_right + m_left) / 2);
      x = -alpha_pade(AmI, m_middle);
      if abs(m_left - m_right) < 2
        m = m_right;
        found = true;
      elseif x < -1 || scalar_error_pade(x, m_middle) >= epsilon
        m_left = m_middle;
      else % (abs(y - log (1 + x)) < epsilon)
        m_right = m_middle;
      end
    end
  end

  function y = alpha(A, k, m)
    %ALPHA    Estimate max(||X^p||^(1/p), ||X^(p+1)||^(1/(p+1))).
    %   ALPHA(A,K,M) estimates max(||X^p||^(1/p), ||X^(p+1)||^(1/(p+1)))
    %   where p is the largest integer such that p(p-1) < K+M+1.
    R = double(A);
    d = floor((1 + sqrt(2 * (m + k) + 5)) / 2);
    dp1 = normAm(R, d)^(1/d);
    d = d + 1;
    dp2 = normAm(R, d)^(1/d);
    y = mp(max(dp1, dp2));
  end

  function S = logm_taylor(A,m)
    %LOGM_TAYLOR   Taylor approximation to matrix logarithm.
    %   LOGM_TAYLOR(A,M) computes the Taylor approximant of degree M to
    %   LOG(EYE(SIZE(A))+A) using Paterson-Stockmeyer algorithm.
    v_0m = cast(0:1:m, class(A));
    c = (1./v_0m').*(-1).^(v_0m' + 1);
    c(1) = zeros(1, class(A));

    ss = ceil(sqrt(m));
    rr = floor(m/ss);

    % Compute first ss+1 powers.
    mpowers = zeros([size(A), ss+1], class(A));
    I = eye(size(A,1), class(A));
    mpowers(:,:,1) = I;
    mpowers(:,:,2) = A;
    for i=3:ss+1
      mpowers(:,:,i) = A * mpowers(:,:,i-1);
    end

    % Evaluate rr-1 polynomials of degree at most ss.
    B = zeros([size(A),ss+1], class(A));
    for kk = 0 : rr-1
      tt = c(ss*kk+1) * I;
      for j=1:ss-1
        tt = tt + c(ss*kk+j+1) * mpowers(:,:,j+1);
      end
      B(:,:,kk+1) = tt;
    end

    % Evaluate last polynomial of degree m mod ss.
    B(:,:,rr+1) = c(m+1) * mpowers(:,:,m-ss*rr+1);
    for j=m-1:-1:ss*rr
      if j == ss*rr
        B(:,:,rr+1) = B(:,:,rr+1) + c(ss*rr+1)*I;
      else
        B(:,:,rr+1) = B(:,:,rr+1) +...
          c(j+1) * mpowers(:,:,m-ss*rr-(m-j)+1);
      end
    end

    % Use Horner scheme to evaluate polynomials.
    As = mpowers(:,:,ss+1);
    S = zeros(size(A,1), class(A));
    for kk=rr:-1:0
      S = S * As + B(:,:,kk+1);
    end
  end

%% Recomputation of diagonal blocks.

  function F = powerm2by2(T, p)
    %POWERM2BY2    Off-diagonal of power of 2-by-2 upper triagnular matrix.
    %   POWERM2BY2(T, P) computes the off-diagonal element of the Pth power
    %   of the 2-by-2 upper triangular matrix T.

    % Equations (5.5) and (5.6) of
    % N. J. Hihgam and L. Lin. A Schur-Pade algorithm for fractional powers
    % of a matrix. SIAM J. Matrix Anal. & App., DOI: 10.1137/10081232X.
    dT = diag(T);
    diag_diff = dT(2) - dT(1);
    dTp = dT.^p;
    logdT = log(dT);
    if diag_diff == 0 % Repeated eigenvalue.
      F = p * T(1,2) * dT(1)^(p-1);
    elseif any(abs(dT) < abs(flip(dT))/2)
      F =  T(1,2) * (dTp(2) - dTp(1)) / diag_diff;
    else                % To avoid numerical cancellation.
      pi_typed = cast(mp('pi'), class(T));
      tmp = p * (atanh(diag_diff / (dT(2) + dT(1))) +...
        1i * pi_typed * unwinding(logdT(2) - logdT(1)));
      F = 2 * T(1,2) * exp(p * (logdT(2) + logdT(1)) / 2) *...
        sinh(tmp) / diag_diff;
    end
  end

  function F = logm2by2(T)
    %LOGM2BY2    Principal logarithm of 2-by-2 upper triangular matrix.
    %   LOGM2BY2(T) computes the logarithm of the 2-by-2 upper triangular
    %   matrix T.

    % Equation (11.28) of
    % N. J. Higham. Functions of Matrices. SIAM, 2008.
    dT = diag(T);
    diag_diff = dT(2) - dT(1);
    logdT = log(dT);
    F = diag(logdT);
    if diag_diff == 0 % Repeated eigenvalue.
      F(1,2) = T(1,2) / dT(1);
    elseif any(abs(dT) < abs(flip(dT))/2)
      F(1,2) =  T(1,2) * (logdT(2) - logdT(1)) / diag_diff;
    else                % To avoid numerical cancellation.
      pi_typed = cast(mp('pi'), class(T));
      tmp = (2 * atanh(diag_diff / (dT(2) + dT(1))) +...
        2 * pi_typed * 1i * unwinding(logdT(2) - logdT(1)));
      F(1,2) = T(1,2) * tmp / diag_diff;
    end
  end

  function r = sqrt_power_1(a,k)
    %SQRT_POWER_1    Accurate computation of a^(1/2^k) - 1.
    %   SQRT_POWER_1(A,N) computes A^(1/2^K)-1 where A is a vector and K is
    %   a natural number.

    % Algorithm 2 from
    % A. H. Al-Mohy.  A more accurate Briggs method for the logarithm.
    % Numer. Algorithms, DOI: 10.1007/s11075-011-9496-z.
    if k == 0
      r = a - 1;
    elseif k > 0
      khat = k;
      if imag(a) < 0
        a = sqrt(a);
        khat = k - 1;
      end
      z0 = a - 1;
      a = sqrt(a);
      r = 1 + a;
      for i = 1 : khat - 1
        a = sqrt(a);
        r = r .* (1 + a);
      end
      r = z0 ./ r;
    else
      error('K must be a natural number.')
    end
  end

  function u = unwinding(z)
    %UNWINDING    Unwinding number.
    %   UNWINDING(Z) computes the unwinding number of the complex number Z.
    pi_typed = cast(mp('pi'), class(z));
    u = ceil((imag(z) - pi_typed) / (2 * pi_typed));
  end

%% Multiprecision Computing Toolbox utility functions.

  function e = myeps(curr_class)
    %COMP_EPS    Machine epsilon.
    %   COMP_EPS(CLASS) computes the machine epsilon of class CLASS.
    if(strcmp(curr_class, 'mp'))
      e = mp('eps');
    else
      e = eps(curr_class);
    end
  end

%% Miscellanea

  function isavailable = check_tb_available(tb_name)
    %CHECK_TB_AVAILABLE    Check availability of toolbox.
    %    CHECK_TB_AVAILABLE(TOOLBOX) checks whether a toolbox whose name
    %    matches exactly TOOLBOX is available.
    v = ver;
    if any(strcmp(tb_name, {v.Name}))
      isavailable = true;
    else
      isavailable = false;
    end
  end

end
