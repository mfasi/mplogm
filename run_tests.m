format compact
paused = false;

% create output directories
mkdir('figs')
mkdir('pngfigs')
mkdir('tabs')

% addpath
addpath('./external')
addpath('./include')
addpath('./matlog')

% input ids
ids_min = 1;
ids_max = 64;

% Sizes
msize = 4;

% Markers
marker_logagm = '^';
marker_logm = 'o';
marker_logmct = 'v';

marker_logp = '*';
marker_logfp = 'x';
marker_logt = 'd';
marker_logft = 's';

% Colours
color_cond =   [0.0 0.9 0.9];
color_logm =   [1.0 0.4 0.4];
color_logmct = [0.0 0.8 0.0];
color_logagm = [0.6 0.0 0.9];

color_logt =   [1.0 0.5 0.0];
color_logft =  [0.0 0.6 0.0];
color_logp =   [0.0 0.0 0.8];
color_logfp =  [0.0 0.0 0.0];

% Line styles
% ls_cond = '-';
% ls_logm = '-';
% ls_logmct = ':';
% ls_logagm = '--';
% ls_logp = ':';
% ls_logt = '-';
% ls_logfp = '-.';
% ls_logft = ':';
ls_cond = '-';
ls_logm = '-';
ls_logmct = '-';
ls_logagm = '-';
ls_logp = '-';
ls_logt = '-';
ls_logfp = '-';
ls_logft = '-';
lw = 1; % linewidth

% precisions
p_vec = [6, 8, 10];
ymins = [1e-66, 1e-258, mp('1e-1026')];
ymaxs = [1e-48, 1e-243, mp('1e-1008')];

%% Figure 3.1
fprintf('*** Figure 3.1: ');
n = 10;
p = 16;
t_global = tic();
test_alpha_bound %%%%%
time = toc(t_global);
fprintf('%.2f s\n', time);

if (paused)
    pause
end

%% Figure 6.1
fprintf('*** Figure 6.1: ');
n = 10;
p = 16;
t_global = tic();
test_bck_fwd_hist %%%%%
time = toc(t_global);
fprintf('%.2f s\n', time);

if (paused)
    pause
end

%% Figure 6.2
fprintf('*** Figure 6.2: ');
t_global = tic();
test_abs_rel_error %%%%%
time = toc(t_global);
fprintf('%.2f s\n', time);

if (paused)
    pause
end

%% Figures 6.3
fprintf('*** Figure 6.3: ');
t_global = tic();
n = 10;
test_logm
time = toc(t_global);
fprintf('%.2f s\n', time);

return

%% Table 6.1
fprintf('*** Table 6.1: ');
t_global = tic();
test_profile_table
time = toc(t_global);
fprintf('%.2f s\n', time);