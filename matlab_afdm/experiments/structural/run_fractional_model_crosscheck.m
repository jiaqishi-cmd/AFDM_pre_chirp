function results = run_fractional_model_crosscheck(options)
%RUN_FRACTIONAL_MODEL_CROSSCHECK Compare two fractional-Doppler channel models.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1
        options = struct();
    end

    N = get_option(options, 'N', 32);
    betaValues = get_option(options, 'beta_values', [0, 0.1, 0.25, 0.4, 0.5]);
    delays = get_option(options, 'delays', [0, 2]);
    integerDopplers = get_option(options, 'integer_dopplers', [0, 2]);
    gains = get_option(options, 'gains', [1, exp(1i * 0.37)] / sqrt(2));
    c1 = get_option(options, 'c1', 7 / (2 * N));
    c2 = get_option(options, 'c2', sqrt(2) / (10 * N));
    samplingFrequency = get_option(options, 'sampling_frequency', N);
    referenceRoot = get_option(options, 'reference_root', default_reference_root(rootDir));

    referenceFunction = fullfile(referenceRoot, 'cconvLTVChannel.m');
    if exist(referenceFunction, 'file') ~= 2
        error(['Reference cconvLTVChannel.m not found. Supply options.reference_root. ' ...
            'Expected default: %s'], referenceRoot);
    end
    addpath(referenceRoot);

    validateattributes(delays, {'numeric'}, {'vector', 'integer', 'nonnegative'});
    if numel(delays) ~= numel(integerDopplers) || numel(delays) ~= numel(gains)
        error('delays, integer_dopplers, and gains must have the same length.');
    end

    transform = build_daft_matrix(N, c1, c2);
    currentVsReference = zeros(numel(betaValues), 1);
    currentVsExplicit = zeros(numel(betaValues), 1);
    referenceVsExplicit = zeros(numel(betaValues), 1);

    for betaIdx = 1:numel(betaValues)
        dopplers = integerDopplers;
        dopplers(end) = dopplers(end) + betaValues(betaIdx);

        [hReference, ~] = cconvLTVChannel( ...
            N, delays / samplingFrequency, ...
            dopplers * samplingFrequency / N, gains, ...
            samplingFrequency, true, c1);
        hExplicit = build_explicit_time_channel(N, c1, gains, delays, dopplers);
        hEffectiveReference = transform * hReference * transform';
        hEffectiveExplicit = transform * hExplicit * transform';
        hEffectiveCurrent = afdm.rx.estimate_effective_channel( ...
            N, c1, c2, gains, delays, dopplers / N);

        currentVsReference(betaIdx) = relative_error( ...
            hEffectiveCurrent, hEffectiveReference);
        currentVsExplicit(betaIdx) = relative_error( ...
            hEffectiveCurrent, hEffectiveExplicit);
        referenceVsExplicit(betaIdx) = relative_error( ...
            hEffectiveReference, hEffectiveExplicit);
    end

    results = table( ...
        betaValues(:), currentVsReference, currentVsExplicit, referenceVsExplicit, ...
        'VariableNames', {'fractional_doppler', 'current_vs_reference', ...
        'current_vs_explicit', 'reference_vs_explicit'});

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    writetable(results, fullfile(outputDir, 'fractional_model_crosscheck.csv'));
    save(fullfile(outputDir, 'fractional_model_crosscheck.mat'), ...
        'results', 'N', 'c1', 'c2', 'delays', 'integerDopplers', 'gains');

    fprintf('Fractional-Doppler model cross-check\n');
    fprintf('Reference: %s\n', referenceFunction);
    disp(results);
end

function transform = build_daft_matrix(N, c1, c2)
    n = (0:N-1).';
    lambda1 = diag(exp(-1i * 2 * pi * c1 * n.^2));
    lambda2 = diag(exp(-1i * 2 * pi * c2 * n.^2));
    transform = lambda2 * (dftmtx(N) / sqrt(N)) * lambda1;
end

function hTime = build_explicit_time_channel(N, c1, gains, delays, dopplers)
    identityN = eye(N);
    n = (0:N-1).';
    hTime = zeros(N, N);

    for pathIdx = 1:numel(gains)
        delay = delays(pathIdx);
        shift = circshift(identityN, delay, 1);
        doppler = diag(exp(-1i * 2 * pi * dopplers(pathIdx) * n / N));
        cppDiagonal = ones(N, 1);
        if delay > 0
            prefixIndices = (0:delay-1).';
            cppDiagonal(1:delay) = exp(-1i * 2 * pi * c1 * ...
                (N^2 - 2 * N * (delay - prefixIndices)));
        end
        hTime = hTime + gains(pathIdx) * diag(cppDiagonal) * doppler * shift;
    end
end

function value = relative_error(first, second)
    value = norm(first - second, 'fro') / max(norm(second, 'fro'), eps);
end

function path = default_reference_root(rootDir)
    afdmParent = fileparts(fileparts(rootDir));
    path = fullfile(afdmParent, 'research_april', 'tmp', ...
        'doubly-dispersive-channel-simulation');
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
