function cfg = config()
%% CONFIG Configuration structure for the Log-SNRAS master pipeline.
%
%   CFG = CONFIG() returns a structure containing all fixed parameters,
%   ephemeris data, file paths, threshold values, and ground-truth keyword
%   lists required to reproduce the Log-SNRAS manuscript.
%
%   The configuration is based on the ephemeris of pi Men c (Huang et al.
%   2018) for TIC 261136679 segments, and on independent literature
%   keywords for all other targets, ensuring zero contamination by the
%   metrics themselves.
%
%   FIELDS
%   ------
%   T0_btjd, Period, T_dur, T_margin
%       Transit ephemeris: T0, period, duration, and time margin for
%       pi Men c (BTJD = BJD_UTC - 2457000).
%   tess_folder, bbs_folder
%       Relative paths to the TESS .fits and BBS .csv/.xlsx input folders.
%   tier1_cut, tier2_cut, threshold_detection, B_boot
%       Thresholds for Tier classification, detection cut-off, and
%       number of bootstrap replicates.
%   artifact_keywords, complex_keywords, other_confirmed
%       Cell arrays of case-insensitive substrings used to identify
%       artifact targets, complex/ambiguous cases, and confirmed planets
%       from independent literature.
%
%   OUTPUT
%   ------
%   cfg : struct
%       Structure with all configuration variables.
%

    % --- 0. EPHEMERIS (pi Men c) ---
    cfg.T0_btjd   = 1425.789204;   % BTJD = BJD_UTC - 2457000
    cfg.Period    = 6.2678399;     % days
    cfg.T_dur     = 0.091;         % days (transit duration)
    cfg.T_margin  = cfg.T_dur / 2 + 0.05;
    
    % --- 1. CONFIGURATION ---
    cfg.tess_folder = 'TESS';
    cfg.bbs_folder  = 'BBS';
    cfg.tier1_cut = 15;
    cfg.tier2_cut = 60;
    cfg.threshold_detection = 7.1;
    cfg.B_boot = 2000;
    
    % --- 2. GROUND-TRUTH KEYWORD LISTS ---
    cfg.artifact_keywords = {'12644769', '5812701', 'v723_mon', 'tic14444029', ...
                              '8462852', '11253226', '3858884', 'proxima', ...
                              '10935310', '6206751', '4544587', '9851944', ...
                              '8912308', '6070714', '9700322', ...
                              '11446443', '5385723', '7943602', '3832716'};
    cfg.complex_keywords  = {'12557548', 'wd_1145', '10001893', '8120608', ...
                              '229747848', '441462308', '278683844', '144440290'};
    cfg.other_confirmed   = {'kepler_1625', 'kepler-1625', 'toi_201', 'toi-201'};
end
