function [A, sA] = getgallery(k, n, par)
% Function test_getall, (probably) written by Bruno Iannazzo
% copied here and just renamed 

    if nargin < 3
        parkms=0.5;
        parpei=0.5;
    else
        parkms=par;
        parpei=par;
    end

    % skipped 'compar','house','integerdata','lauchli','neumann',
    %         'normaldata','uniformdata','wathen','wilk'
    switch k

      case 1, sA='cauchy';A=gallery(sA,n); 
      case 2, sA='chebspec';A=gallery(sA,n); 
      case 3, sA='chebvand';A=gallery(sA,n); 
      case 4, sA='circul';A=gallery(sA,n); 
      case 5, sA='clement';A=gallery(sA,n); 
      case 6, sA='condex';A=gallery(sA,n); 
      case 7, sA='dorr';A=full(gallery(sA,n)); 
      case 8, sA='dramadah';A=gallery(sA,n); 
      case 9, sA='fiedler';A=gallery(sA,n); 
      case 10, sA='forsythe';A=gallery(sA,n); 

      case 11, sA='frank';A=gallery(sA,n); 
      case 12, sA='gcdmat';A=gallery(sA,n); 
      case 13, sA='grcar';A=gallery(sA,n); 
      case 14, sA='hanowa';A=gallery(sA,n); 
      case 15, sA='invhess';A=gallery(sA,n); 
      case 16, sA='kahan';A=gallery(sA,n); 
      case 17, sA='kms';A=gallery(sA,n); 
      case 18, sA='lesp';A=gallery(sA,n); 
      case 19, sA='lotkin';A=gallery(sA,n); 
      case 20, sA='minij';A=gallery(sA,n); 

      case 21, sA='moler';A=gallery(sA,n); 
      case 22, sA='parter';A=gallery(sA,n); 
      case 23, sA='pei';A=gallery(sA,n);
      case 24, sA='prolate';A=gallery(sA,n); 
      case 25, sA='qmult';A=gallery(sA,n);
      case 26, sA='randcolu';A=gallery(sA,n);
      case 27, sA='randcorr';A=gallery(sA,n); 
      case 28, sA='randhess';A=gallery(sA,n); 
      case 29, sA='randjorth';A=gallery(sA,n); 
      case 30, sA='rando';A=gallery(sA,n); 

      case 31, sA='randsvd';A=gallery(sA,n); 
      case 32, sA='riemann';A=gallery(sA,n);
      case 33, sA='ris';A=gallery(sA,n); 
      case 34, sA='sampling';A=gallery(sA,n);
      case 35, sA='toeppd';A=gallery(sA,n); 
      case 36, sA='toeppen';A=full(gallery(sA,n)); 
      case 37, sA='tridiag';A=full(gallery(sA,n)); 
      case 38, A=magic(n);sA='magic';
      case 39, A=hilb(n);sA='hilb';
      case 40, A=pascal(n);sA='pascal';

      case 41, A=wilkinson(n);sA='wilkinson';

        % symbolic lambert not well computed
      case 42, sA='invol';A=gallery(sA,n); % repeated eigenvalues
      case 43, sA='binomial';A=gallery(sA,n); % repeated eigenvalues
      case 44, sA='redheff';A=gallery(sA,n);  % repeated eigenvalues   
      case 45, sA='smoke';A=gallery(sA,n); % 1 eigenvalue is not computed accurately
      case 46, A=invhilb(n);sA='invhilb';

        % singular
      case 47, sA='chow';A=gallery(sA,n); % non diag
      case 48, sA='cycol';A=gallery(sA,n);
      case 49, sA='gearmat';A=gallery(sA,n); 

        % non diagonalizable or similar
      case 50, sA='ipjfact';A=gallery(sA,n); 
      case 51, sA='jordbloc';A=gallery(sA,n); 
      case 52, sA='orthog';A=gallery(sA,n); 
      case 53, sA='triw';A=gallery(sA,n,0.1); 


      case 100, theta=0:pi/8:2*pi-pi/8;
        % eigenvalues around the circle |z-1/2|=3/2 
        v=[2+1e-14*randn(1,14) 1 4];
        M=randn(16);
        A=M*diag(v)/M;
        sA='Eigs on |z-1/2|=3/2';
      case 101, theta=0:pi/8:2*pi-pi/8;
        % eigenvalues around the circle |z+1/3|=1/2 
        v=[-1/6+1e-14*randn(1,14) 1 4];
        M=randn(16);
        A=M*diag(v)/M;
        sA='Eigs on |z+1/2|=1/3';
      case 102,
        D=diag([1 -1/exp(1)+0.1i -1/exp(1)-0.2i 1i/2 2-1i]);
        A=D;
        sA='Diagonal matrix';
      case 103,
        M=randn(5);
        D=diag([1 -1/exp(1)+0.1i -1/exp(1)-0.2i 1i/2 2-1i]); 
        A=M\D*M;
        sA='Random matrix?';
      case 104,
        d=1e-2;D=diag([1 2 -1+1i*d -1-1i*d]);M=rand(4);A=M\D*M;
        sA='Eigs near to the branch cut';
    end

    % s=sprintf('\n\nProcessing now: %s',sA);
    % disp(s);

