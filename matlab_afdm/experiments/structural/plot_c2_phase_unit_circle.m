% PLOT_C2_PHASE_UNIT_CIRCLE
% Visualize quadratic pre-chirp phases on the unit circle for baseline, GPS,
% and proposed c2 constructions.

rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

N = 64;
V = 4;
c2_base = sqrt(2) / (10 * N);
delta_ratio = 0.2;
delta = delta_ratio * c2_base;
gps_pattern = [2 2 1 1];
proposed_pattern = [2 2 1 1];

m = (0:N-1).';
c2_baseline = c2_base * ones(N, 1);
[c2_gps, ~] = afdm.chirp.build_gps_pattern(N, V, gps_pattern);
[c2_proposed, ~] = afdm.chirp.build_proposed_pattern(N, V, proposed_pattern, c2_base, delta);

methods = struct( ...
    'name', {'Traditional AFDM', 'GPS', 'Proposed'}, ...
    'c2', {c2_baseline, c2_gps, c2_proposed}, ...
    'color', {[0.16 0.32 0.70], [0.82 0.28 0.20], [0.12 0.55 0.35]});

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

fig = figure('Name', 'c2 quadratic phase unit circle', 'Color', 'w', ...
    'Position', [100 100 1200 420]);
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

theta = linspace(0, 2*pi, 512);
for idx = 1:numel(methods)
    nexttile;
    z = exp(-1i * 2 * pi * methods(idx).c2(:) .* (m.^2));
    phases = mod(angle(z), 2*pi);
    uniqueBins = numel(unique(round(phases, 10)));

    plot(cos(theta), sin(theta), 'Color', [0.72 0.72 0.72], 'LineWidth', 1.2);
    hold on;
    scatter(real(z), imag(z), 42, methods(idx).color, 'filled', ...
        'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.25);
    axis equal;
    xlim([-1.15 1.15]);
    ylim([-1.15 1.15]);
    grid on;
    xlabel('Real');
    ylabel('Imag');
    title(sprintf('%s\\nunique phase bins: %d / %d', ...
        methods(idx).name, uniqueBins, N));
end

sgtitle(sprintf('Quadratic phase exp(-j2\\pi c_{2,m}m^2), N=%d, V=%d, proposed \\delta/c_2=%.2g', ...
    N, V, delta_ratio), 'FontWeight', 'bold');

pngPath = fullfile(outputDir, ['c2_phase_unit_circle_' timestamp '.png']);
figPath = fullfile(outputDir, ['c2_phase_unit_circle_' timestamp '.fig']);
saveas(fig, pngPath);
savefig(fig, figPath);
fprintf('Saved c2 phase unit circle figure to %s\n', pngPath);
