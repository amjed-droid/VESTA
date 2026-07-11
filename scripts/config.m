function cfg = config()
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
