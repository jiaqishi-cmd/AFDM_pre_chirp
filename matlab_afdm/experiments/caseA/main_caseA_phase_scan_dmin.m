% MAIN_CASEA_PHASE_SCAN_DMIN
% Fixed Case A phase scan for GPS-specific minimum-distance vulnerability.
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

N = 64;
V = 4;
M = N / V;
alpha_max = 3;
c1 = (2 * alpha_max + 1) / (2 * N);
c2_base = sqrt(2) / (10 * N);
theta_list = linspace(0, 2 * pi, 720);

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

gps_pattern = [1 2 1 2];
bestPatternPath = fullfile(outputDir, 'focused_bemani_gps_best_cases.mat');
if exist(bestPatternPath, 'file')
    loaded = load(bestPatternPath, 'focused_best_cases');
    if isfield(loaded, 'focused_best_cases') && isfield(loaded.focused_best_cases, 'mismatch_case')
        gps_pattern = parse_pattern(loaded.focused_best_cases.mismatch_case.gps_pattern);
    end
end

fprintf('========== Case A phase scan dmin ==========\n');
fprintf('N=%d V=%d M=%d c1=%.10g c2_base=%.10g\n', N, V, M, c1, c2_base);
fprintf('Case A: path1=(0,0), path2=(2,2), L=16\n');
fprintf('GPS pattern: [%s]\n', sprintf('%d ', gps_pattern));

deltaOptions.num_random_bpsk = 2000;
deltaOptions.num_random_qpsk = 2000;
deltaOptions.seed = 20260508;
delta_set = generate_delta_set_for_dmin(N, M, deltaOptions);
fprintf('Delta candidates: %d\n', size(delta_set.delta, 2));

scanOptions.N = N;
scanOptions.V = V;
scanOptions.alpha_max = alpha_max;
scanOptions.c1 = c1;
scanOptions.c2_base = c2_base;
scanOptions.gps_pattern = gps_pattern;
scanOptions.proposed_pattern = gps_pattern;
scan = scan_theta_dmin_caseA(theta_list, delta_set, scanOptions);

phase_scan_table = scan.phase_scan_table;
best_case = scan.best_case; %#ok<NASGU>
best_theta = best_case.theta; %#ok<NASGU>
best_pattern = gps_pattern; %#ok<NASGU>
best_delta_type_GPS = best_case.delta_type_GPS; %#ok<NASGU>
best_dmin2_base = best_case.dmin2_base; %#ok<NASGU>
best_dmin2_GPS = best_case.dmin2_GPS; %#ok<NASGU>
best_dmin2_prop = best_case.dmin2_prop; %#ok<NASGU>
best_delta_GPS = scan.best_delta_GPS; %#ok<NASGU>

save(fullfile(outputDir, 'caseA_phase_scan_dmin_results.mat'), ...
    'phase_scan_table', 'best_case', 'best_theta', 'best_pattern', ...
    'best_delta_type_GPS', 'best_dmin2_base', 'best_dmin2_GPS', ...
    'best_dmin2_prop', 'best_delta_GPS', 'scan');
writetable(phase_scan_table, fullfile(outputDir, 'caseA_phase_scan_dmin_results.csv'));

print_phase_scan_summary(phase_scan_table);
fprintf('\nSelected best theta: %.8f rad = %.6f*pi\n', best_case.theta, best_case.theta_over_pi);
fprintf('dmin2 base %.4g | GPS %.4g | prop %.4g | gaps %.2f/%.2f dB | delta %s\n', ...
    best_case.dmin2_base, best_case.dmin2_GPS, best_case.dmin2_prop, ...
    best_case.gap_base_GPS_dB, best_case.gap_prop_GPS_dB, best_case.delta_type_GPS);

plot_caseA_phase_scan_results(scan, outputDir);
fprintf('Saved Case A phase scan results to %s\n', outputDir);

function print_phase_scan_summary(T)
    print_sorted(T, 'dmin2_GPS', 'ascend', 'Top 20 by smallest GPS dmin2');
    print_sorted(T, 'ratio_GPS_base', 'ascend', 'Top 20 by smallest GPS/base ratio');
    print_sorted(T, 'gap_base_GPS_dB', 'descend', 'Top 20 by largest baseline-GPS gap');
    valid = T.dmin2_base > 0.1 & T.dmin2_prop > 0.1;
    if ~any(valid)
        valid = T.dmin2_base > 0.03 & T.dmin2_prop > 0.03;
    end
    if any(valid)
        print_sorted(T(valid, :), 'score', 'descend', 'Top 20 with baseline/proposed not too small');
    else
        print_sorted(T, 'score', 'descend', 'Top 20 by score');
    end
end

function print_sorted(T, fieldName, direction, titleText)
    [~, order] = sort(T.(fieldName), direction);
    top = T(order(1:min(20, height(T))), :);
    fprintf('\n%s:\n', titleText);
    disp(top(:, {'theta','theta_over_pi','dmin2_base','dmin2_GPS','dmin2_prop', ...
        'ratio_GPS_base','gap_base_GPS_dB','gap_prop_GPS_dB', ...
        'delta_type_GPS','score'}));
end

function pattern = parse_pattern(text)
    if isnumeric(text)
        pattern = text;
        return;
    end
    switch char(text)
        case 'pattern_alt1'
            pattern = [1 2 1 2];
        case 'pattern_alt2'
            pattern = [2 1 2 1];
        case 'pattern_half1'
            pattern = [1 1 2 2];
        case 'pattern_half2'
            pattern = [2 2 1 1];
        case 'pattern_all1'
            pattern = [1 1 1 1];
        otherwise
            nums = regexp(char(text), '\d+', 'match');
            pattern = str2double(nums);
    end
end
