function results = validate_gps_cycle_closed_form(options)
%VALIDATE_GPS_CYCLE_CLOSED_FORM Compare closed form with numerical cycles.

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
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    delayList = get_option(options, 'delay_list', 1:8);
    dopplerList = get_option(options, 'doppler_list', -alphaMax:alphaMax);

    modulationCases = struct( ...
        'name', {"BPSK", "QPSK", "16QAM"}, ...
        'M_mod', {2, 4, 16}, ...
        'mod_type', {'psk', 'psk', 'qam'}, ...
        'rotational_order', {2, 4, 4});

    c2 = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'rational');
    h1 = afdm.analysis.build_path_matrix(N, c1, c2, 0, 0);
    numRows = numel(modulationCases) * numel(delayList) * numel(dopplerList);
    rows = repmat(empty_row(), numRows, 1);
    cursor = 0;

    for modulationIdx = 1:numel(modulationCases)
        modulation = modulationCases(modulationIdx);
        deltaSet = afdm.delta.generate_sparse_difference_set( ...
            1, modulation.M_mod, modulation.mod_type, 1);
        for delay2 = delayList
            for doppler2 = dopplerList
                predicted = afdm.analysis.gps_cycle_risk_condition( ...
                    N, c1, delay2, doppler2, modulation.rotational_order);
                h2 = afdm.analysis.build_path_matrix( ...
                    N, c1, c2, delay2, doppler2);
                numerical = afdm.analysis.monomial_cycle_analysis( ...
                    h1' * h2, deltaSet.difference_alphabet, 1e-12);
                numericalRisk = numerical.num_compatible_eigenvectors > 0;
                exactPhaseCondition = ...
                    afdm.analysis.monomial_rotational_cycle_condition( ...
                    h1' * h2, modulation.rotational_order, 1e-10);

                cursor = cursor + 1;
                row = empty_row();
                row.modulation = modulation.name;
                row.delay2 = delay2;
                row.doppler2 = doppler2;
                row.eta = predicted.eta;
                row.cycle_count = predicted.cycle_count;
                row.cycle_length = predicted.cycle_length;
                row.phase_closure_residue = predicted.phase_closure_residue;
                row.simple_congruence_risk = ...
                    predicted.predicts_compatible_cycle;
                row.predicted_risk = ...
                    exactPhaseCondition.predicts_compatible_cycle;
                row.numerical_risk = numericalRisk;
                row.match = row.predicted_risk == row.numerical_risk;
                rows(cursor) = row;
            end
        end
    end

    summary = struct2table(rows);
    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'pattern', pattern, 'delay_list', delayList, ...
        'doppler_list', dopplerList);
    results.summary = summary;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'gps_cycle_closed_form_validation.mat'), 'results');
    writetable(summary, ...
        fullfile(outputDir, 'gps_cycle_closed_form_validation.csv'));

    fprintf('Closed-form matches: %d/%d\n', sum(summary.match), height(summary));
    for modulationIdx = 1:numel(modulationCases)
        selected = summary.modulation == modulationCases(modulationIdx).name;
        fprintf('%-6s predicted=%d numerical=%d mismatches=%d\n', ...
            char(modulationCases(modulationIdx).name), ...
            sum(summary.predicted_risk(selected)), ...
            sum(summary.numerical_risk(selected)), ...
            sum(~summary.match(selected)));
    end
    if ~all(summary.match)
        disp(summary(~summary.match, :));
        error('Closed-form condition did not match numerical cycle analysis.');
    end
end

function row = empty_row()
    row = struct( ...
        'modulation', "", ...
        'delay2', NaN, ...
        'doppler2', NaN, ...
        'eta', NaN, ...
        'cycle_count', NaN, ...
        'cycle_length', NaN, ...
        'phase_closure_residue', NaN, ...
        'simple_congruence_risk', false, ...
        'predicted_risk', false, ...
        'numerical_risk', false, ...
        'match', false);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
