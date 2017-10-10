function [Aexpm, set, id, nmats] = logm_testmats(k, n)

    rng(0);

    % Bad: 1, 17, 23,
    expm_ids = [4, 6, 8, 9, 10, 11, 12, 13, 14,... % 17, 23, 25
        15, 16, 18, 19, 21, 22, 24, 26, 27,...
        28, 30]; % 22
    expm_nmats = length(expm_ids); % removed 51
    gallery1_ids = [1, 3, 8, 10, 13, 14, 16, 17, 18, 19,...
        22, 24, 25, 26, 27, 28, 31, 33, 36, 37,...
        39, 41, 44, 45, 47, 48, 49, 52, 53]; % 29 (51)
    gallery1_nmats = length(gallery1_ids);
    gallery2_ids = [6, 7, 12, 13, 15, 17, 20, 22, 23, 27, 35, 11]; % 11 (62)
    gallery2_nmats = length(gallery2_ids);
    mymats = 3;
    nmats = expm_nmats + gallery1_nmats + gallery2_nmats + mymats;

    if (k <= expm_nmats)
        set = 1;
        id = expm_ids(k);
        A = expm_testmats(id, n);
        thelog = true;
    elseif (k <= expm_nmats + gallery1_nmats)
        set = 2;
        id = gallery1_ids(k-expm_nmats);
        A = getgallery(id, n);
        thelog = true;
    elseif (k <= expm_nmats + gallery1_nmats + gallery2_nmats)
        set = 3;
        id = gallery2_ids(k-expm_nmats-gallery1_nmats);
        A = getgallery(id, n);
        thelog = false;
    elseif (k <= nmats)
        set = 4;
        offset = expm_nmats + gallery1_nmats + gallery2_nmats;
        switch k
            case offset + 1.
                id = 101;
                A = [0.01, 0.95; 0, 0.04];

            case offset + 2.
                % Referee report.
                id = 102;
                alpha = 0.1;
                beta = 1e2;
                A = exp(alpha) * [1 beta; 0 1];

            case offset + 3.
               % Referee report.
               id = 103;
               a = 0.1;
               b = 1e3;
               c = 1e3;
               A = exp(a) * [1 b b^2/2 + c; 0 1 b; 0 0 1];
        end 
        thelog = false;
    else
        error('This gallery contains only %d matrices.', nmats);
    end

    if (thelog)
        Aexpm = expm(A);
    else
        Aexpm = A;
    end

end
