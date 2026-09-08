function results = run_gps_path_support_grid(options)
%RUN_GPS_PATH_SUPPORT_GRID Scan constellation-compatible path cycles.

    if nargin < 1
        options = struct();
    end

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    N = get_option(options, 'N', 64);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 3);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    delta = get_option(options, 'proposed_delta', c2Base / 16);
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    delayList = get_option(options, 'delay_list', 1:3);
    dopplerList = get_option(options, 'doppler_list', -alphaMax:alphaMax);
    differenceAlphabet = get_option(options, ...
        'difference_alphabet', [-2i; 2i]);
    resultLabel = get_option(options, 'label', 'bpsk');

    schemes = ["baseline", "GPS", "proposed"];
    c2Values = { ...
        c2Base, ...
        afdm.chirp.build_gps_pattern(N, V, pattern), ...
        afdm.chirp.build_proposed_pattern(N, V, pattern, c2Base, delta)};
    h1 = cell(numel(schemes), 1);
    for schemeIdx = 1:numel(schemes)
        h1{schemeIdx} = afdm.analysis.build_path_matrix( ...
            N, c1, c2Values{schemeIdx}, 0, 0);
    end

    numRows = numel(delayList) * numel(dopplerList) * numel(schemes);
    rows = repmat(empty_row(), numRows, 1);
    cursor = 0;

    for delay2 = delayList
        for doppler2 = dopplerList
            for schemeIdx = 1:numel(schemes)
                h2 = afdm.analysis.build_path_matrix( ...
                    N, c1, c2Values{schemeIdx}, delay2, doppler2);
                relativeOperator = h1{schemeIdx}' * h2;
                cycleAnalysis = afdm.analysis.monomial_cycle_analysis( ...
                    relativeOperator, differenceAlphabet);

                cursor = cursor + 1;
                row = empty_row();
                row.delay2 = delay2;
                row.doppler2 = doppler2;
                row.scheme = schemes(schemeIdx);
                row.effective_shift = cycleAnalysis.constant_shift;
                row.gcd_shift_N = gcd(round(cycleAnalysis.constant_shift), N);
                row.cycle_length = N / row.gcd_shift_N;
                row.compatible_eigenvectors = ...
                    cycleAnalysis.num_compatible_eigenvectors;
                row.min_compatible_weight = ...
                    cycleAnalysis.min_compatible_weight;
                rows(cursor) = row;
            end
        end
    end

    summary = struct2table(rows);
    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'delta', delta, 'pattern', pattern, ...
        'delay_list', delayList, 'doppler_list', dopplerList);
    results.summary = summary;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputTag = ['gps_path_support_grid_' sanitize_label(resultLabel)];
    save(fullfile(outputDir, [outputTag '.mat']), 'results');
    writetable(summary, fullfile(outputDir, [outputTag '.csv']));
    plot_grid(results, outputDir, outputTag);
    print_summary(summary);
end

function plot_grid(results, outputDir, outputTag)
    gpsRows = results.summary(results.summary.scheme == "GPS", :);
    delayList = results.config.delay_list;
    dopplerList = results.config.doppler_list;
    gridValues = NaN(numel(delayList), numel(dopplerList));
    for rowIdx = 1:height(gpsRows)
        delayIdx = find(delayList == gpsRows.delay2(rowIdx), 1);
        dopplerIdx = find(dopplerList == gpsRows.doppler2(rowIdx), 1);
        gridValues(delayIdx, dopplerIdx) = ...
            gpsRows.min_compatible_weight(rowIdx);
    end

    figure('Name', 'GPS path-support cycle risk', 'Color', 'w');
    imagesc(dopplerList, delayList, gridValues);
    set(gca, 'YDir', 'normal');
    colorbar;
    xlabel('Second-path integer Doppler');
    ylabel('Second-path delay');
    title('Minimum BPSK-compatible cycle weight for GPS');
    saveas(gcf, fullfile(outputDir, [outputTag '.png']));
end

function label = sanitize_label(value)
    label = regexprep(char(string(value)), '[^A-Za-z0-9_-]', '_');
end

function print_summary(T)
    schemes = unique(T.scheme, 'stable');
    for schemeIdx = 1:numel(schemes)
        subset = T(T.scheme == schemes(schemeIdx), :);
        risky = subset.compatible_eigenvectors > 0;
        fprintf('%-9s risky supports %d/%d', ...
            char(schemes(schemeIdx)), sum(risky), height(subset));
        if any(risky)
            fprintf(', compatible weights [%s]', ...
                sprintf('%d ', unique(subset.min_compatible_weight(risky))));
        end
        fprintf('\n');
    end
end

function row = empty_row()
    row = struct( ...
        'delay2', NaN, ...
        'doppler2', NaN, ...
        'scheme', "", ...
        'effective_shift', NaN, ...
        'gcd_shift_N', NaN, ...
        'cycle_length', NaN, ...
        'compatible_eigenvectors', NaN, ...
        'min_compatible_weight', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
