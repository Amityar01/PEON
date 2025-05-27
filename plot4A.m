
%% plot4A.m
% Analyzes and plots the population-level firing rates of identified PEONs
% (Probability Encoding Omission Neurons) in response to omissions, the 
% Omission Preferred (OP) tone, and the Omission Non-Preferred (ONP) tone
% across different OP tone probability conditions. This script generates 
% figures corresponding to Figure 4A of the manuscript (Yaron et al., 2025),
% showing both response curves and Spearman correlation coefficients.
%
% Workflow:
% 1. Takes PEON indices ('PEONs_training') identified using split-half analysis
%    (training on odd trials) and their preferred tone directions 
%    ('preferred_tone_direction') as input, along with 3D matrices of
%    trial-by-trial responses ('allommat', 'allAmat', 'allBmat') and
%    'testing_indices' (even trials) for independent validation.
% 2. For each PEON and each of the 8 OP tone probability conditions (0%-95%):
%    a. Calculates mean omission response from 'allommat' using 'testing_indices',
%       aligning data based on 'preferred_tone_direction' to ensure OP tone
%       omissions show positive correlations with probability.
%    b. Calculates mean response to the OP tone (from 'allAmat' or 'allBmat'
%       depending on preferred direction) for the 7 conditions where tones are present.
%    c. Calculates mean response to the ONP tone similarly, using appropriate
%       probability indexing based on the PEON's preferred direction.
% 3. Aggregates individual PEON responses by probability condition across the
%    population to calculate means and standard errors for plotting.
% 4. Performs Spearman correlations between individual neuron responses and
%    their corresponding OP tone probabilities for each response type.
% 5. Generates two publication-ready figures:
%    a. Line plot showing population mean firing rates (±SEM) vs OP tone
%       probability for omission, OP tone, and ONP tone responses (Fig 4A top).
%    b. Bar plot showing Spearman correlation coefficients for the three
%       response types with significance markers (Fig 4A bottom).
% 6. Exports data to Excel files for reproducibility and further analysis.
%
% Key Findings Illustrated:
%   - PEONs show selective omission responses that increase with OP tone probability
%   - Tone responses are broad and inversely related to probability (adaptation)
%   - Asymmetry between selective prediction errors and broad sensory responses
%
% Inputs (expected in MATLAB workspace, typically from 'FindPEONS.m' & 'allmats.m'):
%   - PEONs_training: Vector of indices for neurons classified as PEONs using
%                     split-half analysis (identified on odd trials).
%   - preferred_tone_direction: Vector (-1 or 1) for each PEON, determining
%                               which physical tone (A or B) is the OP tone.
%   - allommat: 3D matrix (neurons x 50_trials x 8_probability_conditions) of
%               trial-by-trial omission responses (from 'allmats.m').
%   - allAmat: 3D matrix (neurons x 50_trials x 7_tone_conditions) for Tone A
%              responses across probability conditions (from 'allmats.m').
%   - allBmat: 3D matrix (neurons x 50_trials x 7_tone_conditions) for Tone B
%              responses across probability conditions (from 'allmats.m').
%   - testing_indices: Vector of trial indices (e.g., 2:2:50 for even trials)
%                      used for independent validation of PEONs identified on odd trials.
%
% Outputs:
%   - Figure 1: Line plot of population mean firing rates (±SEM) vs OP tone
%               probability for omission, OP tone, and ONP tone responses.
%               Saved as 'PEON_responses_analysis.pdf'.
%   - Figure 2: Bar plot of Spearman correlation coefficients with significance
%               markers for the three response types. Shows the asymmetric
%               encoding: positive correlation for omissions, negative for OP tones,
%               positive for ONP tones (reflecting deviance detection).
%   - 'Figure_4A_data.xlsx': Excel file containing:
%     * Response_Data sheet: Population means and SEMs for each probability condition
%     * Correlations sheet: Spearman correlation statistics
%     * Individual_Points sheet: All individual neuron responses for reproducibility
%   - Console output: Spearman correlation results (rho, p-values, sample sizes)
%
% Dependencies:
%   - Standard MATLAB statistics and plotting functions
%   - Data structures from 'FindPEONS.m' and 'allmats.m'
%   - Excel writing capability for data export
%
% Note: This script demonstrates the key asymmetry in PEON responses:
%       selective encoding of negative prediction errors (omissions) while
%       maintaining broad sensory responses to both tones, supporting the
%       lateral prediction suppression model proposed in the manuscript.
%
% Author: Amit Yaron
% Corresponding to: Figure 4A in Yaron et al. (2025)

peon_idx = PEONs_training;
nPeons = length(peon_idx);
probability_values = [0, 0.05, 0.1, 0.2, 0.75, 0.85, 0.90, 0.95];

% Initialize arrays for collecting responses separately for each probability
om_by_prob = cell(1,8);
op_by_prob = cell(1,8);
onp_by_prob = cell(1,8); 

% For each neuron
for i = 1:nPeons
    neuron = peon_idx(i);
    
    % For each probability condition
    for p = 1:8
        % For omissions:
        if preferred_tone_direction(i) == 1
            om_resp = mean(allommat(neuron, testing_indices, p));
        else
            om_resp = mean(allommat(neuron, testing_indices, 9-p));
        end
        om_by_prob{p} = [om_by_prob{p}; om_resp];
        
        % For tones - we need to handle boundary cases
        if p >= 2 && p <= 8  % For OP tone (conditions 2-8)
            if preferred_tone_direction(i) == 1
                % If preferred direction is 1, OP is tone B
                op_resp = mean(allBmat(neuron, testing_indices, p-1));
            else
                % If preferred direction is -1, OP is tone A
                op_resp = mean(allAmat(neuron, testing_indices, 9-p));
            end
            op_by_prob{p} = [op_by_prob{p}; op_resp];
        end
        
        if p >= 1 && p <= 7  % For ONP tone (conditions 1-7)
            if preferred_tone_direction(i) == 1
                % If preferred direction is 1, ONP is tone A
                onp_resp = mean(allAmat(neuron, testing_indices, p));
            else
                % If preferred direction is -1, ONP is tone B
                onp_resp = mean(allBmat(neuron, testing_indices, 8-p));
            end
            onp_by_prob{p} = [onp_by_prob{p}; onp_resp];
        end
    end
end

% Calculate means and errors
means_om = zeros(1,8);
err_om = zeros(1,8);
means_op = NaN(1,8);  % Use NaN for positions without data
err_op = NaN(1,8);
means_onp = NaN(1,8);
err_onp = NaN(1,8);

for p = 1:8
    % Omissions - all positions
    means_om(p) = mean(om_by_prob{p});
    err_om(p) = std(om_by_prob{p})/sqrt(length(om_by_prob{p}));
    
    % OP - positions 2-8
    if p >= 2 && p <= 8
        means_op(p) = mean(op_by_prob{p});
        err_op(p) = std(op_by_prob{p})/sqrt(length(op_by_prob{p}));
    end
    
    % ONP - positions 1-7
    if p >= 1 && p <= 7
        means_onp(p) = mean(onp_by_prob{p});
        err_onp(p) = std(onp_by_prob{p})/sqrt(length(onp_by_prob{p}));
    end
end

% For correlation analysis - flatten data
all_om_resp = [];
all_om_prob = [];
all_op_resp = [];
all_op_prob = [];
all_onp_resp = [];
all_onp_prob = [];

for p = 1:8
    % Omissions
    all_om_resp = [all_om_resp; om_by_prob{p}];
    all_om_prob = [all_om_prob; repmat(probability_values(p), length(om_by_prob{p}), 1)];
    
    % OP tone
    if p >= 2 && p <= 8
        all_op_resp = [all_op_resp; op_by_prob{p}];
        all_op_prob = [all_op_prob; repmat(probability_values(p), length(op_by_prob{p}), 1)];
    end
    
    % ONP tone
    if p >= 1 && p <= 7
        all_onp_resp = [all_onp_resp; onp_by_prob{p}];
        all_onp_prob = [all_onp_prob; repmat(probability_values(p), length(onp_by_prob{p}), 1)];
    end
end

% Calculate correlations
[rho_om, p_om] = corr(all_om_prob, all_om_resp, 'Type', 'Spearman');
[rho_op, p_op] = corr(all_op_prob, all_op_resp, 'Type', 'Spearman');
[rho_onp, p_onp] = corr(all_onp_prob, all_onp_resp, 'Type', 'Spearman');

% Display results
fprintf('Omission responses: r = %.4f, p = %.2e, n = %d\n', rho_om, p_om, length(all_om_resp));
fprintf('Preferred tone (OP) responses: r = %.4f, p = %.2e, n = %d\n', rho_op, p_op, length(all_op_resp));
fprintf('Non-preferred tone (ONP) responses: r = %.4f, p = %.2e, n = %d\n', rho_onp, p_onp, length(all_onp_resp));

% Create figure for plotting
figure('Position', [100, 100, 600, 800]);

% Top subplot - mean responses with error bars
subplot(2,1,1);
hold on;

errorbar(1:8, means_om, err_om, 'o-', 'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'DisplayName', 'Omission response');
errorbar(1:8, means_op, err_op, 's-', 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0 0.4470 0.7410], 'DisplayName', 'Tone response to O_P');
errorbar(1:8, means_onp, err_onp, 'd-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Tone response to O_{NP}');

xlabel('Probability of O_P tone in sequence (%)', 'FontSize', 12);
ylabel('Firing rate (spk/s)', 'FontSize', 12);
xticks(1:8);
xticklabels({'0', '5', '10', '20', '75', '85', '90', '95'});
grid on;
legend('Location', 'best');
hold off;

% Bottom subplot - correlation bar plot
figure
bar_values = [rho_op, rho_onp, rho_om];
bar_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0 0 0];
b = bar(1:3, bar_values, 'FaceColor', 'flat');
for k = 1:3
    b.CData(k,:) = bar_colors(k,:);
end
hold on;

% Add significance markers
for i = 1:3
    p_val = [p_op, p_onp, p_om];
    if p_val(i) < 0.05
        text(i, bar_values(i) + 0.05, '*', 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
end

line([0.5 3.5], [0 0], 'Color', 'k', 'LineStyle', '--');
ylim([-.5 .5]);
xticks(1:3);
xticklabels({'Response to O_P', 'Response to O_{NP}', 'Omission response'});
ylabel('Spearman rho', 'FontSize', 12);
grid on;

% Export figure
set(gcf, 'Renderer', 'painters');
exportgraphics(gcf, 'PEON_responses_analysis.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

%%
% Extract Figure 4A data
prob_values = [0, 5, 10, 20, 75, 85, 90, 95];

% Create data table for the three response types
figure_4A_data = table();
figure_4A_data.Probability_Percent = prob_values';
figure_4A_data.Omission_Mean = means_omission';
figure_4A_data.Omission_SEM = err_omission';
figure_4A_data.TonePref_Mean = means_tonePref';
figure_4A_data.TonePref_SEM = err_tonePref';
figure_4A_data.ToneNonPref_Mean = means_toneNonPref';
figure_4A_data.ToneNonPref_SEM = err_toneNonPref';

% Save main data
writetable(figure_4A_data, 'Figure_4A_data.xlsx', 'Sheet', 'Response_Data');

% Save correlation statistics
stats_4A = {'Omission_rho', rho_omission; 'Omission_p', p_omission; ...
           'Preferred_rho', rho_pref; 'Preferred_p', p_pref; ...
           'NonPreferred_rho', rho_nonpref; 'NonPreferred_p', p_nonpref};
writetable(cell2table(stats_4A, 'VariableNames', {'Statistic', 'Value'}), ...
           'Figure_4A_data.xlsx', 'Sheet', 'Correlations');

fprintf('Figure 4A data saved to Figure_4A_data.xlsx\n');
%%
%% --- Build long-format table of individual responses for Figure 4A ---

% ---- prerequisites -----------------------------------------------
%  allom_peon_mean :  nPeons x 8   (omission, cols = [0 5 10 20 75 85 90 95])
%  allA_peon_mean  :  nPeons x 7   (preferred-tone,   cols = [5 10 20 75 85 90 95])
%  allB_peon_mean  :  nPeons x 7   (non-pref tone,    cols = [0 5 10 20 75 85 90])
%  neuronIDs       :  nPeons x 1 cell array, e.g. {'N001','N002',…}
% -------------------------------------------------------------------

prob_all = [0 5 10 20 75 85 90 95];   % 8 probabilities (col order in omission)

rows = {};    % will become a cell array of rows
% ------------------------------------------------------------------
% Make simple text labels like N0001, N0002, … in the PEON order
% ------------------------------------------------------------------
neuronIDs = arrayfun(@(k) sprintf('N%04d', k), peon_idx, ...
                     'UniformOutput', false);   % 123×1 cell array

for i = 1:nPeons
    % ----- Omission (8 values) ------------------------------------
    for j = 1:8
        rows(end+1,:) = { neuronIDs{i}, prob_all(j), 'Omission', ...
                          allom_peon_mean(i,j) };
    end

    % ----- Preferred tone  (7 values: prob index 2–8) --------------
    for j = 2:8
        rows(end+1,:) = { neuronIDs{i}, prob_all(j), 'TonePref', ...
                          allA_peon_mean(i,j-1) };
    end

    % ----- Non-preferred tone (7 values: prob index 1–7) ----------
    for j = 1:7
        rows(end+1,:) = { neuronIDs{i}, prob_all(j), 'ToneNonPref', ...
                          allB_peon_mean(i,j) };
    end
end

fig4A_indiv = cell2table(rows, ...
        'VariableNames', {'NeuronID','Probability_Percent','ResponseType','FiringRate'});

% Append to the same workbook
writetable(fig4A_indiv, 'Figure_4A_data.xlsx', ...
           'Sheet', 'Individual_Points');

fprintf('Added %d individual rows to sheet "Individual_Points".\n',height(fig4A_indiv));
