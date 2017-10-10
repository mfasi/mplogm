function [varargout] = logm_mp_full(A, varargin)
%LOGM_MP_FULL  Transformation-free multiprecision matrix logarithm.
%   This code requires the Advanpix Multiprecision Computing Toolbox
%   (see www.advanpix.com).
%
%   [X, S, M] = logm_mp_full(A) computes the principal matrix logarithm X
%   of the matrix A. The output parameters S and M are the number of square
%   roots taken by the algorithm and the degree of the Pade approximation,
%   respectively.
%
%   [...] = logm_mp_full(...,'precision',DIGITS) specifies the number of
%   digits to be used in the computation. Default is mp.Digits() if A
%   is of class mp. The computation is performed in single or double
%   arithmetic if A is of class single or double, respectively, and
%   the precision is not specified.
%
%   [...] = logm_mp_full(...,'epsilon',EPSILON) specifies the tolerance to
%   be used to evaluate the Pade approximations. Default is machine epsilon
%   of the precision of A if A is of class 'single' or 'double', mp.eps()
%   if A is of class 'mp'.
%
%   [...] = logm_mp_full(...,'maxsqrtm', MAXSQRTM) specifies the maximum
%   number of square roots to be taken. Default is 100.
%
%   [...] = logm_mp_full(...,'maxdegree', MAXDEGREE) specifies the maximum
%   degree of the Pade approximants. Default is 100.
%
%   [...] = logm_mp_full(...,'approx', KIND), where KIND='diagonal' or
%   KIND='taylor', specifies the kind of Pade approximant to use. Default
%   is diagonal. The efficient evaluation of truncated Taylor series
%   requires storing several additional matrices having the same size as A.
%
%   Reference: Algorithm 5.1 in
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
    'Multiprecision Computing Toolbox,\nwhich does not appear'...
    'to be installed on this system.\nSee www.advanpix.com.'];
  error('logm_mp:missing_toolbox', err_msg);
end

%% Parse and validate input.
p = inputParser;
addRequired(p, 'A', @(x)(ismatrix(x) && isnumeric(x)));
addParameter(p, 'precision', [], @(x)(x == round(x) && x > 0))
addParameter(p, 'epsilon', [], @(x)(x > 0 && x < 1));
addParameter(p, 'maxsqrt', 100, @(x)(x == round(x) && x > 0));
addParameter(p, 'maxdegree', 100, @(x)(x == round(x) && x > 0));
addParameter(p, 'approx', 'diagonal',...
  @(x)(ischar(x) && strcmp(x, 'taylor') || strcmp(x, 'diagonal')));
parse(p,A,varargin{:});

A = p.Results.A;
digits = p.Results.precision;
epsilon = p.Results.epsilon;
maxsqrt = p.Results.maxsqrt;
maxdegree = p.Results.maxdegree;

if strcmp(p.Results.approx, 'taylor')
  usetaylor = true;
  logm_pade = @logm_taylor;
  eval_pade_error = @scalar_error_taylor;
  alpha_pade = @(A, m)(alpha(A,m,0));
else
  usetaylor = false;
  logm_pade = @logm_diagonal;
  eval_pade_error = @scalar_error_diagonal;
  alpha_pade = @(A, m)(alpha(A,m,m));
end
min_pade_deg = @(A, epsilon, maxd)...
  (find_min_deg(A, epsilon, maxd, eval_pade_error));

[n1, n] = size(A);
if n ~= n1
  error('The matrix ``A'' must be square.');
end

%% Parse and validate output.
if nargout > 3
  error('This function returns at most three values.');
end

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

A0 = A;
e = eig(A);
use_ccc = false;
if any(imag(e) == 0 & real(e) <= 0)
  use_ccc = true;
  warning(['Principal logarithm not defined for matrices with\n',...
    'nonpositive real eigenvalues, the solution conforms\n',...
    'to the rule of counter-clockwise continuity.']);
end

%% Compute initial square roots.
% Square until norm(T - eye(n)) < maxnorm.
s = 0;
Z = cast([], class(A));
P = cast([], class(A));
id_n = eye(n, class(A));
AmI = A - id_n;
eps_normAmI = epsilon * norm(AmI, 1);
a = -alpha_pade(AmI, maxdegree);
while (norm(AmI, 1) >= 1 ||...
    eval_pade_error(a, maxdegree) >= eps_normAmI &&...
    s < maxsqrt)
  [A, Z, P, s] = compute_sqrtm(A, Z, P, s);
  AmI = A - id_n;
  eps_normAmI = epsilon * norm(AmI, 1);
  a = -alpha_pade(AmI, maxdegree);
end

%% Estimate the initial value for m.
m = min_pade_deg(A, eps_normAmI, maxdegree);

try_sqrtm = true;
while (try_sqrtm && s < maxsqrt)
  if usetaylor
    newm = max(ceil((sqrt(m) - 4)^2)-1, 1);
  else
    newm = max(m - 7, 1);
  end
  a = -alpha_pade(AmI, m) / 2;
  boundabs2 = eval_pade_error(a, newm);

  % Check if squaring and taking lower degree approximant will work.
  if boundabs2 > eps_normAmI       % no
    try_sqrtm = false;
  else                             % yes
    [A, Z, P, s] = compute_sqrtm(A, Z, P, s);
    AmI = A - id_n;
    eps_normAm1 = epsilon * norm(AmI, 1);
    m = min_pade_deg(A, eps_normAm1, m);
  end
end

if s < 2
  Y = A - id_n;
else
  Y = P \ Z;
end

%% Compute the matrix logarithm.
X = logm_pade(Y, m);

% Scale back
X = 2^s * X;
if isreal(A) && ~use_ccc
  X = real(X);
end

%% Prepare output.
[varargout{1}] = X;
k = 2;
if k < nargout
  [varargout{k}] = s;
  k = k + 1;
  [varargout{k}] = m;
end

%% Restore mp working precision.
mp.Digits(current_digits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             SUBFUNCTIONS                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Square root.

  function [X, Z, P, s] = compute_sqrtm(A, Z, P, s)
    %SQRTM_DBP  Matrix square by Denman-Beavers iteration in product form.

    % This function is based on the funnction SQRTM_DBP from
    % Nicholas J. Higham, The Matrix Function Toolbox,
    % http://www.ma.man.ac.uk/~higham/mftoolbox

    nn = length(A);
    I = eye(nn);

    tol = sqrt(length(A)) * myeps('mp')/2;
    X = A;
    M = A;
    it = 1;
    maxit = 50;

    scale = 1;
    while(norm(M - I,'fro') > tol && it <= maxit)

      if scale == 1
        g = (abs(det(M)))^(-1/(2*nn));
        X = g*X;
        M = g^2 * M;
      end

      Xold = X;
      invM = inv(M);

      X = X*(I + invM)/2;
      M = 0.5*(I + (M + invM)/2);

      reldiff = norm(X - Xold, 'fro') / norm(X, 'fro');
      if reldiff < 1e-2
        scale = 0;
      end  % Switch to no scaling.
      it = it + 1;

    end
    if it == maxit
      warning('Not converged after %2.0f iterations', maxit)
    end
    s = s + 1;
    switch s
      case 1
        Z = X - I;
      case 2
        P = I + X;
      otherwise
        P = P * (I + X);
    end
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

  function m = find_min_deg(A, epsilon, maxd, scalar_error_pade)
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
      x = -alpha(AmI, m_middle, m_middle);

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

  function S = logm_diagonal(A, m)
    %LOGM_PADE   Pade approximation to matrix logarithm.
    %   LOGM_PADE(A,M) computes the [M/M] Pade approximant to
    %   LOG(EYE(SIZE(A))+A) using the partial fraction expansion.
    [nodes, wts] = mp.GaussLegendre(m, 0, 1);
    nodes = cast(nodes, class(A));
    wts = cast(wts, class(A));
    S = zeros(size(A, 1), class(A));
    for j = 1:m
      S = S + wts(j) * ((eye(size(A, 1)) + nodes(j) * A) \ A);
    end
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
