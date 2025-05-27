%% Neuronal Response Analysis to Probabilistic Stimuli - PEON Identification
% This script analyzes neuronal responses to auditory stimuli presented with
% varying probabilities to identify Probability Encoding Omission Neurons (PEONs).
% It uses a split-half approach (ODD vs. EVEN trials) for robustness and
% performs statistical tests to determine significance.
%
% Data dimensions assumed for input 'data' (e.g., allommat):
%   n_neurons x n_trials_per_prob x n_probabilities
%
% Variable names reverted closer to original for easier comparison.

%% Section 1: Setup and Data Preparation ==================================
% --- User-Defined Parameters ---
ANALYSIS_MODE = 'ODD'; % Options: 'ODD', 'EVEN', 'ALL'. Determines which trials are used.
ALPHA_CORR = 0.05;     % Significance level for Spearman correlation
ALPHA_WILCOX = 0.05;   % Significance level for Wilcoxon signed-rank test

% --- Load and Prepare Data ---
deep=T.depth(1:8:end)';
clear is_response_significantal
omind =find(deep<1800 & deep>0);
allmats
data = allommat; % Using omission responses for PEON identification

% Get depths (example, ensure T and relevant fields exist)
% deep = T.depth(1:8:end)'; % Original line - check if correct for your T table
% omind = find(deep < 1800 & deep > 0); % Example filtering by depth

% --- Define Data Dimensions and Constants ---
[num_neurons, num_trials_per_prob, num_probabilities] = size(data);
TOTAL_TRIALS = num_trials_per_prob; % Should be 50 based on original code
PROBABILITY_CONDITIONS = 1:num_probabilities; % Index for conditions (1 to 8)
PROBABILITY_VALUES = [0, 5, 10, 20, 75, 85, 90, 95]; % Actual probability percentages
HIGH_PROB_CONDITIONS_IDX = 5:8; % Indices for high probability conditions (corresponding to 75% to 95% in aligned data)

fprintf('--- Starting PEON Analysis ---\n');
fprintf('Analysis Mode: %s\n', ANALYSIS_MODE);
fprintf('Total Neurons: %d\n', num_neurons);
fprintf('Trials per Probability: %d\n', TOTAL_TRIALS);
fprintf('Number of Probabilities: %d\n', num_probabilities);
fprintf('Alpha (Correlation): %.2f\n', ALPHA_CORR);
fprintf('Alpha (Wilcoxon): %.2f\n', ALPHA_WILCOX);
fprintf('----------------------------------\n');

%% Section 2: Initial Analysis - All Data (For Figure 2A) ================
% Calculates Spearman correlation for all neurons using the entire dataset.
% Used for the histogram in Figure 2A showing distribution of rhos.

fprintf('Section 2: Calculating correlations on ALL data...\n');
rho_all = zeros(num_neurons, 1);       % Absolute rho values (magnitude) - Reverted Name
p_values_all = zeros(num_neurons, 1);  % p-values - Reverted Name
rho_alldir = zeros(num_neurons, 1);    % Sign of rho (+1 or -1) - Reverted Name

% Prepare data for correlation: flatten trials and probabilities
prob_labels_all = repmat(PROBABILITY_CONDITIONS, 1, TOTAL_TRIALS); % Repeat prob index (1-8) for each trial
% Reshape data: num_neurons x (num_probabilities * num_trials_per_prob)
data_flat_all = reshape(permute(data, [1, 3, 2]), num_neurons, []);

% Calculate correlations for each neuron
for neuron_idx = 1:num_neurons
    [rho_temp, p_values_all(neuron_idx)] = corr(prob_labels_all', data_flat_all(neuron_idx, :)', ...
        'Type', 'Spearman', 'Rows', 'complete'); % Handle potential NaNs
    if ~isnan(rho_temp) % Handle potential NaN results from corr
        rho_all(neuron_idx) = abs(rho_temp); % Store magnitude
        rho_alldir(neuron_idx) = sign(rho_temp); % Store direction
    else
        rho_all(neuron_idx) = NaN;
        rho_alldir(neuron_idx) = NaN;
        p_values_all(neuron_idx) = NaN;
    end
end
fprintf('Section 2: Done.\n');

%% Section 3: Split-Half Analysis Setup ==================================
% Defines trial indices based on ANALYSIS_MODE and prepares datasets for
% training (PEON identification) and testing (plotting/validation).

fprintf('Section 3: Setting up analysis for mode: %s...\n', ANALYSIS_MODE);

switch ANALYSIS_MODE
    case 'ODD'
        % Use ODD trials for training, EVEN for testing
        training_indices = 1:2:TOTAL_TRIALS; % e.g., [1, 3, ..., 49]
        testing_indices = 2:2:TOTAL_TRIALS;  % e.g., [2, 4, ..., 50]
        fprintf('Using ODD trials for Training, EVEN trials for Testing.\n');
    case 'EVEN'
        % Use EVEN trials for training, ODD for testing
        training_indices = 2:2:TOTAL_TRIALS;
        testing_indices = 1:2:TOTAL_TRIALS;
        fprintf('Using EVEN trials for Training, ODD trials for Testing.\n');
    case 'ALL'
        % Use ALL trials for training, no separate testing set defined here
        training_indices = 1:TOTAL_TRIALS;
        testing_indices = []; % No separate test set needed for PEON ID steps
        fprintf('Using ALL trials for Training.\n');
    otherwise
        error('Invalid ANALYSIS_MODE specified. Choose ODD, EVEN, or ALL.');
end

num_training_trials = length(training_indices);
num_testing_trials = length(testing_indices); % Will be 0 if ANALYSIS_MODE is 'ALL'

% Create training dataset (used for identifying PEONs)
training_data = data(:, training_indices, :);

% Create testing dataset (used for plotting final PEONs, unless mode is 'ALL')
if ~isempty(testing_indices)
    testing_data = data(:, testing_indices, :);
else
    testing_data = []; % Explicitly empty if mode is 'ALL'
end

% Prepare flattened training data for correlation analysis
training_flat = reshape(permute(training_data, [1, 3, 2]), num_neurons, []);

% Create probability labels corresponding to the training set
prob_labels_train = repmat(PROBABILITY_CONDITIONS, 1, num_training_trials);

fprintf('Section 3: Done.\n');

%% Section 4: Identify Initial PEON Candidates (Training Set) ============
% Calculates correlations on the training set data to find neurons showing
% significant modulation by probability.

fprintf('Section 4: Identifying initial PEON candidates using Training data...\n');
rho_training = NaN(num_neurons, 1); % Use NaN for initialization
p_values_training = NaN(num_neurons, 1);

% Calculate Spearman correlation for each neuron using training data
for neuron_idx = 1:num_neurons
    [rho_training(neuron_idx), p_values_training(neuron_idx)] = ...
        corr(prob_labels_train', training_flat(neuron_idx, :)', 'Type', 'Spearman', 'Rows', 'complete');
end

% Identify initial candidates based on correlation significance
% NOTE: Using PEONs_training for the *initial* list here, will be overwritten in Sec 6
PEONs_training_initial_idx = find(p_values_training < ALPHA_CORR); % find automatically handles NaNs
num_initial_candidates = length(PEONs_training_initial_idx);

fprintf('Found %d initial PEON candidates (p < %.2f on Training data).\n', num_initial_candidates, ALPHA_CORR);
fprintf('Section 4: Done.\n');

%% Section 5: Response Analysis and Alignment ============================
% Aligns responses for all neurons based on the preferred direction of
% correlation (sign of rho_training from Section 4). This ensures that for
% all neurons, probability condition 1 represents the lowest probability
% of the preferred omission, and condition 8 represents the highest.

fprintf('Section 5: Aligning responses based on preferred direction (from Training rho)...\n');

% Use rho_training from Section 4 to determine preferred direction for *all* neurons
preferred_tone_direction = sign(rho_training); % Reverted Name
% Handle cases where rho might be exactly zero or NaN (assign default direction, e.g., +1)
preferred_tone_direction(isnan(preferred_tone_direction) | preferred_tone_direction == 0) = 1;

% --- Align TRAINING Data (needed for Wilcoxon test in Section 6) ---
% Preallocate aligned training data matrix for all neurons
% Using PEON_responses_orderedtr name for closer match to original concept
PEON_responses_orderedtr = NaN(size(training_data)); % Reverted Name (conceptually)

for neuron_idx = 1:num_neurons
    neuron_training_responses = squeeze(training_data(neuron_idx, :, :)); % num_training_trials x num_probabilities

    if preferred_tone_direction(neuron_idx) == -1
        % Flip probabilities dimension if preferred direction is negative
        PEON_responses_orderedtr(neuron_idx, :, :) = fliplr(neuron_training_responses);
    else
        % Keep original order if preferred direction is positive
        PEON_responses_orderedtr(neuron_idx, :, :) = neuron_training_responses;
    end
end

% --- Align TESTING Data (needed for plotting/analysis later) ---
% Preallocate aligned testing data matrix for all neurons
if ~isempty(testing_data)
    % Using PEON_responses_ordered name for closer match to original concept
    PEON_responses_ordered = NaN(size(testing_data)); % Reverted Name (conceptually)
    for neuron_idx = 1:num_neurons
        neuron_testing_responses = squeeze(testing_data(neuron_idx, :, :)); % num_testing_trials x num_probabilities

        if preferred_tone_direction(neuron_idx) == -1
            PEON_responses_ordered(neuron_idx, :, :) = fliplr(neuron_testing_responses);
        else
            PEON_responses_ordered(neuron_idx, :, :) = neuron_testing_responses;
        end
    end
else
    % If mode is 'ALL', testing data will be used from training data later if needed
    PEON_responses_ordered = [];
end

fprintf('Section 5: Done.\n');


%% Section 6: Refined PEON Identification (Wilcoxon Test) ================
% Refines the PEON list by adding a criterion: the neuron must show a
% significant response (Wilcoxon signed-rank test vs 0, two-tailed) in the
% high probability conditions of the *aligned training data*.
% **MODIFIED:** Reverted signrank call to default two-tailed test.

fprintf('Section 6: Refining PEON list with Wilcoxon test on aligned Training data...\n');

% Perform Wilcoxon test for *all* neurons using aligned training data
wilcox_p_values_all = NaN(num_neurons, 1); % Still store p-value for reference if needed
% Using is_response_significantal name for closer match to original concept
is_response_significantal = false(num_neurons, 1); % Reverted Name (conceptually)

for neuron_idx = 1:num_neurons
    % Get aligned training responses for high probability conditions
    % Use PEON_responses_orderedtr (aligned training data)
    high_prob_responses = PEON_responses_orderedtr(neuron_idx, :, HIGH_PROB_CONDITIONS_IDX);
    high_prob_responses_flat = high_prob_responses(:); % Flatten to 1D array

    % Remove NaNs before testing
    high_prob_responses_flat = high_prob_responses_flat(~isnan(high_prob_responses_flat));

    % Check if data is suitable for signrank (not empty, not all zeros)
    if ~isempty(high_prob_responses_flat) && any(high_prob_responses_flat ~= 0)
        % Perform Wilcoxon signed-rank test against zero (TWO-TAILED by default)
        % Get p-value; h output is not used here.
        [p_wilcox, ~] = signrank(high_prob_responses_flat, 0); % Use default two-tailed test
        wilcox_p_values_all(neuron_idx) = p_wilcox; % Store p-value
        % Store significance result by comparing p-value to ALPHA_WILCOX
        is_response_significantal(neuron_idx) = (p_wilcox < ALPHA_WILCOX); % Compare p-value
    else
        % Cannot perform test if data is empty, all zeros, or contains NaNs only
        wilcox_p_values_all(neuron_idx) = NaN;
        is_response_significantal(neuron_idx) = false;
    end
end

% Define final PEONs: Must meet BOTH criteria
% 1. Significant correlation on training data (p < ALPHA_CORR from Section 4)
% 2. Significant response in high-prob conditions (Wilcoxon p < ALPHA_WILCOX from this section)
% NOTE: Overwriting PEONs_training with the final list, matching original script's likely behaviour
PEONs_training = find( (p_values_training < ALPHA_CORR) & is_response_significantal ); % Reverted Name, uses p-value comparison via is_response_significantal
num_final_peons = length(PEONs_training); % Use final PEON list length

% Identify other categories for pie chart
correlated_not_significant_idx = find( (p_values_training < ALPHA_CORR) & ~is_response_significantal );
% Neurons where correlation p-value was >= alpha OR where p-value was NaN
% **MODIFIED:** Corrected logic for counting categories.
no_correlation_idx = find( p_values_training >= ALPHA_CORR | isnan(p_values_training) );

num_corr_not_sig = length(correlated_not_significant_idx);
num_no_corr = length(no_correlation_idx);

fprintf('Refined PEON Analysis Results:\n');
fprintf('Initial candidates (Corr p < %.2f): %d\n', ALPHA_CORR, num_initial_candidates);
fprintf('Neurons with Sig. Wilcoxon response (p < %.2f, two-tailed): %d\n', ALPHA_WILCOX, sum(is_response_significantal & ~isnan(wilcox_p_values_all))); % Clarified test type
fprintf('Final PEONs (Corr + Wilcoxon): %d\n', num_final_peons);
fprintf('Correlated but Not Sig. Response: %d\n', num_corr_not_sig);
fprintf('No Correlation: %d\n', num_no_corr); % Changed label slightly
% Sanity check count
if (num_final_peons + num_corr_not_sig + num_no_corr) ~= num_neurons
    warning('Neuron category counts do not sum to total neurons! Check logic. Sum=%d, Total=%d', ...
            num_final_peons + num_corr_not_sig + num_no_corr, num_neurons);
else
     fprintf('Category counts sum correctly to %d.\n', num_neurons);
end
fprintf('Section 6: Done.\n');

% --- Prepare data from final PEONs for plotting ---
% Using mean_responses_validated name for closer match to original concept
if ~isempty(PEONs_training) % Check if final PEON list is not empty
    if ~isempty(PEON_responses_ordered) % Use testing data if available (ODD or EVEN mode)
        % Mean responses across trials for FINAL PEONs using aligned TESTING data
        mean_responses_validated = squeeze(mean(PEON_responses_ordered(PEONs_training, :, :), 2)); % Reverted Name
        data_source_for_plots = 'Testing Data';
    else % Use training data if mode is 'ALL'
        % Mean responses across trials for FINAL PEONs using aligned TRAINING data
        mean_responses_validated = squeeze(mean(PEON_responses_orderedtr(PEONs_training, :, :), 2)); % Reverted Name
        data_source_for_plots = 'Training Data (ALL Mode)';
        fprintf('Warning: Using training data for subsequent plots as ANALYSIS_MODE is ALL.\n');
    end
    % Get preferred direction for final PEONs only
    preferred_direction_final_peons = preferred_tone_direction(PEONs_training); % Use final PEON list
else
    % Handle case where no PEONs are found
    mean_responses_validated = [];
    preferred_direction_final_peons = [];
    data_source_for_plots = 'N/A';
    fprintf('Warning: No final PEONs identified.\n');
end


%% Section 7: Summary Statistics & Plotting ==============================
% Generates summary plots based on the FINAL PEON list.

fprintf('Section 7: Generating summary plots using %s...\n', data_source_for_plots);

% --- Pie Chart ---
figure('Name', 'Neuron Distribution Pie Chart');
pie_data = [num_no_corr, num_corr_not_sig, num_final_peons]; % Use counts calculated in Sec 6
pie_labels = {'No Correlation', 'Correlated Not Sig.', 'PEON'};
% Create labels with numbers
pie_numLabels = cell(1, length(pie_data));
for i = 1:length(pie_data)
    % Avoid displaying label if count is zero
    if pie_data(i) > 0
        pie_numLabels{i} = sprintf('%s: %d', pie_labels{i}, pie_data(i));
    else
        pie_numLabels{i} = ''; % Empty label if count is zero
    end
end
% Define custom colors
pie_colors = [0.2 0.6 0.8;  % Blue
              0.8 0.3 0.5;  % Red
              0.3 0.8 0.3]; % Green
colormap(pie_colors);
% Only plot non-zero segments
non_zero_idx = find(pie_data > 0);
if ~isempty(non_zero_idx)
    pie(pie_data(non_zero_idx), pie_numLabels(non_zero_idx));
else
    text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center'); % Placeholder if all counts are zero
end
title(sprintf('Distribution of Neuron Classification Results (%s)', ANALYSIS_MODE));
set(gcf, 'Renderer', 'painters'); % Improve vector export
exportgraphics(gcf, sprintf('piechart_neuron_distribution_%s.pdf', ANALYSIS_MODE), 'ContentType', 'vector', 'BackgroundColor', 'none');
fprintf('Saved pie chart.\n');

% --- Boxplot of Mean Omission Responses ---
% This plot shows the distribution of the *mean* response per neuron.
% Commented out as requested, replaced by the direct correlation plot in Section 8.
%{
if ~isempty(mean_responses_validated) % Use final PEON mean responses
    figure('Name', 'PEON Mean Omission Response Boxplot', 'Position', [100, 100, 400, 600]);
    % Flatten data for boxplot: num_final_peons * num_probabilities
    aligned_responses_flat_validated = mean_responses_validated(:); % Reverted Name
    % Create grouping variable based on probability condition index (1-8)
    aligned_probabilities_idx_flat = repmat(PROBABILITY_CONDITIONS, num_final_peons, 1);
    aligned_probabilities_idx_flat = aligned_probabilities_idx_flat(:);

    boxplot(aligned_responses_flat_validated, aligned_probabilities_idx_flat, ...
        'Labels', arrayfun(@(x) sprintf('%g%%', x), PROBABILITY_VALUES, 'UniformOutput', false), ...
        'Colors', 'b', 'Notch', 'off', 'Symbol', ''); % Use actual percentages as labels
    hold on;
    % Add means as distinct markers (red diamonds)
    mean_across_final_peons = mean(mean_responses_validated, 1);
    scatter(PROBABILITY_CONDITIONS, mean_across_final_peons, 100, 'r', 'diamond', 'filled');
    plot(PROBABILITY_CONDITIONS, mean_across_final_peons, 'r--', 'LineWidth', 1.5);
    hold off;

    xlabel('Probability of Preferred Tone (%)', 'FontSize', 12);
    ylabel('Mean Omission Response (Aligned, Spikes/s)', 'FontSize', 12);
    title(sprintf('PEON Omission Response vs Probability (%s)', ANALYSIS_MODE));
    line(xlim, [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Zero line
    grid on;
    set(gca, 'FontSize', 11);
    set(gcf, 'Renderer', 'painters');
    exportgraphics(gcf, sprintf('PEONs_mean_response_boxplot_%s.pdf', ANALYSIS_MODE), 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('Saved mean response boxplot.\n');
else
    fprintf('Skipping mean response boxplot (no final PEONs found).\n');
end
%}
fprintf('Section 7: Done (Boxplot of mean responses commented out).\n');


%% Section 8: Direct Correlation Analysis & Final Plot (Using Mean Responses) ===================
% Calculates the Spearman correlation using the *mean response per neuron*
% for the final PEONs (using testing data if available, training data otherwise),
% after alignment based on preferred direction.
% Generates the final boxplot showing mean per neuron, overall means, and stats,
% matching the style requested by the user (image_65c058.png, n=984).
% **MODIFIED:** Removed CI calculation and display from plot text box.

fprintf('Section 8: Calculating direct correlation on MEAN responses and generating final plot...\n');

if ~isempty(mean_responses_validated) % Use mean responses calculated in Sec 6
    data_source_corr = data_source_for_plots; % Source determined in Sec 6

    % --- Prepare Data (Using Means) ---
    % Flatten the mean data: (num_final_peons * num_probabilities) x 1
    all_om_resp_means = mean_responses_validated(:); % Use means

    % Create corresponding probability labels (indices 1-8 first)
    prob_rep1 = repmat(PROBABILITY_CONDITIONS, num_final_peons, 1); % num_final_peons x num_probabilities
    all_prob_idx_flat_direct = prob_rep1(:); % Flatten indices to match responses

    % Map condition indices (1-8) to actual probability values (%) for correlation
    all_om_prob = PROBABILITY_VALUES(all_prob_idx_flat_direct)'; % Use actual values 0-95

    % --- Calculate Direct Spearman Correlation (on Means) ---
    [rho_om, p_om] = corr(all_om_prob, all_om_resp_means, ... % Use means
        'Type', 'Spearman', 'Rows', 'complete'); % Handle NaNs
    n_datapoints = sum(~isnan(all_om_resp_means) & ~isnan(all_om_prob)); % Count valid pairs (should be num_final_peons * 8)

    fprintf('Direct correlation on MEAN responses from %s (Spearman):\n', data_source_corr);
    fprintf('rho = %.4f, p = %.2e, n = %d data points (neurons*conditions)\n', rho_om, p_om, n_datapoints);

    % --- Confidence Interval Calculation Removed ---
    % fprintf('Confidence interval calculation removed as requested.\n');

    % --- Generate Final Plot (Matching image_65c058.png Style - Using Means) ---
    figure('Name', 'Final Omission Response Boxplot with Scatter (Mean Responses)', 'Position', [100, 100, 500, 600]);

    % Create the boxplot using mean responses grouped by probability index (1-8)
    boxplot(all_om_resp_means, all_prob_idx_flat_direct, ...
        'Labels', arrayfun(@(x) sprintf('%g', x), PROBABILITY_VALUES, 'UniformOutput', false), ... % Use values 0-95 for labels
        'Colors', 'b', ...
        'Symbol', '', ... % Hide outlier markers
        'Widths', 0.6);
    hold on;

    % Add individual MEAN data points with jitter
    for p_idx = 1:num_probabilities
        % Get mean data points for this probability condition index
        points = all_om_resp_means(all_prob_idx_flat_direct == p_idx);
        points = points(~isnan(points)); % Remove NaNs for plotting
        if isempty(points); continue; end

        % Create horizontal jitter relative to the index p_idx
        jitter_strength = 0.3; % Adjust jitter amount if needed
        jitter = (rand(size(points)) - 0.5) * jitter_strength;

        % Plot points (mean per neuron) with transparency
        scatter(p_idx + jitter, points, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.4); % Semi-transparent black dots
    end

    % Calculate and plot overall mean across neurons for each probability
    mean_overall = mean(mean_responses_validated, 1); % Mean of means
    scatter(PROBABILITY_CONDITIONS, mean_overall, 100, 'r', 'diamond', 'filled', 'MarkerEdgeColor', 'k'); % Red diamond with black edge
    plot(PROBABILITY_CONDITIONS, mean_overall, 'r--', 'LineWidth', 1.5); % Dashed red line connecting means

    % Add statistics text box (matching image_65c058.png, WITHOUT CI)
    stats_text = sprintf('Spearman \\rho = %.2f\n(p = %.2e, n = %d)', ... % Removed CI line
        rho_om, p_om, n_datapoints);
    % Position the text box (adjust coordinates based on expected y-range)
    ax_limits_y = ylim; % Get current limits after plotting means
    ax_limits_x = xlim;
    text_x_pos = ax_limits_x(1) + 0.65 * diff(ax_limits_x); % Position towards right
    text_y_pos = ax_limits_y(1) + 0.1 * diff(ax_limits_y);  % Position towards bottom
    text(text_x_pos, text_y_pos, stats_text, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'BackgroundColor', [0.9 0.9 0.9], ...
        'EdgeColor', 'k', ...
        'FontSize', 10);

    hold off;
    % Labels and title
    xlabel('Probability of O_P Tone (%)', 'FontSize', 12);
    ylabel('Omission Response (spikes/s)', 'FontSize', 12);
    title('Omission Response by Probability', 'FontSize', 14);

    % Add reference line at y=0
    line(xlim, [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '-', 'LineWidth', 1);

    % Enhance appearance
    set(gca, 'LineWidth', 1, 'FontSize', 10, 'XTickLabelRotation', 45); % Rotate x-labels like image 1
    box on;
    grid off; % Grid was off in image 1

    % Set y-limits (adjust to match image 1 style, e.g., -20 to 40)
    ylim([-25 45]); % Adjust as needed

    % Export figure
    set(gcf, 'Renderer', 'painters');
    exportgraphics(gcf, sprintf('final_omission_MEAN_response_correlation_%s.pdf', ANALYSIS_MODE), 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('Saved final MEAN response correlation boxplot.\n');

else
    fprintf('Skipping final plot (no final PEONs found).\n');
    rho_om = NaN; % Define output vars even if skipped
    p_om = NaN;
    n_datapoints = 0;
    % ci_lower_rho = NaN; % CI removed
    % ci_upper_rho = NaN;
end
fprintf('Section 8: Done.\n');


fprintf('--- PEON Analysis Complete ---\n');

% --- Available Output Variables ---
% PEONs_training: Indices of the final PEONs relative to the original neuron list.
% mean_responses_validated: Mean response across trials for each final PEON for each probability (aligned).
% preferred_direction_final_peons: Preferred direction (+1 or -1) for final PEONs.
% rho_om, p_om: Spearman correlation results from direct analysis on MEAN responses (Section 8).
% n_datapoints: Number of valid data points used for direct correlation (neurons * conditions).
% rho_all, rho_alldir: Used for Figure 2A histogram (Section 2).
% Counts for pie chart: num_final_peons, num_corr_not_sig, num_no_corr (Section 6 results).
% ANALYSIS_MODE: The mode used for this run ('ODD', 'EVEN', 'ALL').

%%
% Then extract the data:
figure_S2A_data = [];
for i = 1:length(all_om_resp_means)
    prob_idx = all_prob_idx_flat_direct(i);
    prob_value = PROBABILITY_VALUES(prob_idx);
    response = all_om_resp_means(i);
    figure_S2A_data = [figure_S2A_data; prob_value, response];
end

% Save S2A data
figure_S2A_table = array2table(figure_S2A_data, 'VariableNames', {'Probability_Percent', 'Individual_Response'});
overall_means_S2A = [PROBABILITY_VALUES', mean_overall'];
overall_means_S2A_table = array2table(overall_means_S2A, 'VariableNames', {'Probability_Percent', 'Overall_Mean'});
stats_S2A = {'Spearman_rho', rho_om; 'p_value', p_om; 'n_datapoints', n_datapoints};
stats_S2A_table = cell2table(stats_S2A, 'VariableNames', {'Statistic', 'Value'});

writetable(figure_S2A_table, 'Figure_S2A_data.xlsx', 'Sheet', 'Individual_Points');
writetable(overall_means_S2A_table, 'Figure_S2A_data.xlsx', 'Sheet', 'Overall_Means');
writetable(stats_S2A_table, 'Figure_S2A_data.xlsx', 'Sheet', 'Statistics');
fprintf('Data saved to Figure_2D_data.xlsx\n');