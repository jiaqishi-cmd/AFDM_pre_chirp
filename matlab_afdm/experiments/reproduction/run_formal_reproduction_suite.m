function run_formal_reproduction_suite()
%RUN_FORMAL_REPRODUCTION_SUITE Recreate the main result figures with formal settings.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    logPath = fullfile(outputDir, ['formal_repro_' timestamp '.log']);
    diary(logPath);
    cleanup = onCleanup(@() diary('off'));

    fprintf('Formal reproduction started at %s\n', datestr(now));
    fprintf('Root: %s\n', rootDir);
    fprintf('Log: %s\n\n', logPath);

    run_step('structural metrics', @() run_script(rootDir, 'experiments/structural/test_c2_structural_metrics.m'));
    run_step('PAPR CCDF comparison, 10000 frames', ...
        @() run_papr_ccdf_comparison(10000, struct()));
    run_step('proposed delta CCDF curves, 10000 frames', ...
        @() run_proposed_delta_ccdf_curves(10000, struct()));
    run_step('proposed delta sweep', ...
        @() run_script(rootDir, 'experiments/papr/run_delta_sweep_proposed.m'));
    run_step('proposed delta BER curves', ...
        @() run_script(rootDir, 'experiments/papr/run_proposed_delta_ber_curves.m'));
    run_step('Choi-style random-channel BER', ...
        @() run_script(rootDir, 'experiments/random_channel/run_choi_style_random_channel_ber.m'));
    run_step('Case A fixed-channel deep BER', ...
        @() run_script(rootDir, 'experiments/caseA/run_caseA_fixed_ber_deep.m'));
    run_step('Case A theta scan', ...
        @() run_script(rootDir, 'experiments/caseA/run_caseA_theta_scan.m'));
    run_step('Case A theta enhanced plots', ...
        @() run_script(rootDir, 'experiments/caseA/plot_caseA_theta_scan_enhanced.m'));

    fprintf('\nFormal reproduction finished at %s\n', datestr(now));
end

function run_script(rootDir, scriptPath)
    cmd = sprintf("clearvars; cd('%s'); setup_paths(pwd); run('%s');", ...
        escape_quotes(rootDir), escape_quotes(scriptPath));
    evalin('base', cmd);
end

function text = escape_quotes(text)
    text = strrep(text, '''', '''''');
end

function run_step(name, fn)
    fprintf('\n===== START %s | %s =====\n', name, datestr(now));
    try
        fn();
        fprintf('===== PASS %s | %s =====\n', name, datestr(now));
    catch ME
        fprintf(2, '===== FAIL %s | %s =====\n', name, datestr(now));
        fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        rethrow(ME);
    end
end
