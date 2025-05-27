%% plotmeanstoneom2.m
% This script analyzes and plots the population-level firing rates of identified
% PEONs (Probability Encoding Omission Neurons) in response to:
%   1. Omissions (of the Omission Preferred - OP - tone).
%   2. Presentation of the OP tone.
%   3. Presentation of the Omission Non-Preferred (ONP) tone.
% Responses are analyzed across 8 OP tone probability conditions. The script
% generates figures corresponding to Figure 4A of the manuscript (Yaron et al.,
% 2025), depicting response curves and Spearman correlations.
%
% Workflow:
% 1. Inputs: PEON indices ('PEONs_training'), their preferred tone directions
%    ('preferred_tone_direction'), 3D trial-by-trial response matrices
%    ('allommat', 'allAmat', 'allBmat'), 'testing_indices' (for split-half),
%    and 'probability_values'.
% 2. For each identified PEON:
%    a. Calculates the mean omission response for each of the 8 OP tone
%       probability conditions using 'allommat' and 'testing_indices', aligning
%       data based on 'preferred_tone_direction'.
%    b. Calculates the mean response to the OP tone for the 7 relevant
%       probability conditions. This involves selecting the correct physical tone
%       (A or B, from 'allAmat'/'allBmat') and appropriate probability indexing
%       based on the PEON's 'preferred_tone_direction'.
%    c. Calculates the mean response to the ONP tone similarly for the 7
%       relevant probability conditions.
% 3. Collates all these individual PEON mean responses.
% 4. Calculates population mean and Standard Error of the Mean (SEM) for omission,
%    OP tone, and ONP tone responses for each of the 8 unique OP tone
%    probability values, for plotting the line graph.
% 5. Performs Spearman correlations between the collected individual PEON mean
%    responses and their corresponding OP tone probabilities for each of the
%    three response types (omission, OP tone, ONP tone).
% 6. Generates and saves two figures:
%    a. A line plot (Displayed): Population mean firing rates (+/- SEM) vs.
%       OP tone probability.
%    b. A bar plot (Saved as 'Fig4A_Correlations.pdf'): Spearman correlation
%       coefficients for the three response types.
%
% Inputs (expected to be in the MATLAB workspace):
%   - PEONs_training: Vector of indices for neurons classified as PEONs
%                     (typically output from 'FindPEONS.m').
%   - preferred_tone_direction: Vector (-1 or 1 for each PEON in PEONs_training)
%                               from 'FindPEONS.m'.
%   - allommat: 3D matrix (neurons x trials x 8 probability_conditions) of
%               trial-by-trial omission responses (from 'allmats.m').
%   - allAmat:  3D matrix (neurons x trials x 7 tone_A_conditions) for Tone A
%               responses (from 'allmats.m').
%   - allBmat:  3D matrix (neurons x trials x 7 tone_B_conditions) for Tone B
%               responses (from 'allmats.m').
%   - testing_indices: Vector of trial indices (e.g., 2:2:50 for EVEN trials)
%                      from 'FindPEONS.m' for calculating mean responses.
%   - probability_values: Vector (1x8) of the actual OP tone probability
%                         percentages (e.g., [0, 0.05, ..., 0.95]).
%
% Outputs:
%   - Figure 1 (Displayed): Line plot of mean responses (Fig 4A top-like).
%   - 'Fig4A_Correlations.pdf' (Saved): Bar plot of correlations (Fig 4A bottom-like).
%   - Console output: Spearman correlation results (rho and p-values).
%   - Workspace variables: means_om_plot, sem_om_plot, etc. (for plotting);
%     rho_om, p_om, rho_pref, p_pref, rho_nonpref, p_nonpref (correlation results).
%
% Dependencies:
%   - Standard MATLAB plotting and statistics functions.
%
% Author: Amit Yaron
peon_idx = PEONs_training;  % PEON indices
nPeons = length(peon_idx);
probability_values = [0, 0.05, 0.1, 0.2, 0.75, 0.85, 0.90, 0.95];

% Initialize arrays for collecting responses and probabilities
all_om_resp = [];
all_om_prob = [];
all_pref_resp = [];
all_pref_prob = [];
all_nonpref_resp = []; 
all_nonpref_prob = [];

% For each neuron
for i = 1:nPeons
    neuron = peon_idx(i);
    
    % Process omission responses for all 8 probability conditions
    for p = 1:8
        if preferred_tone_direction(i) == 1
            % Direction 1: original order
            omission_resp = mean(allommat(neuron, testing_indices, p));
            prob_val = probability_values(p);
        else
            % Direction -1: reversed order
            omission_resp = mean(allommat(neuron, testing_indices, 9-p));
            prob_val = probability_values(p);  % Keep original probability order
        end
        
        all_om_resp = [all_om_resp; omission_resp];
        all_om_prob = [all_om_prob; prob_val];
    end
    
    % Process preferred tone responses
    for p = 1:7
        if preferred_tone_direction(i) == 1
            % Direction 1: preferred tone is B
            pref_resp = mean(allBmat(neuron, testing_indices, p));
            prob_val = probability_values(p+1);  % B appears in conditions 2-8
        else
            % Direction -1: preferred tone is A
            pref_resp = mean(allAmat(neuron, testing_indices, p));
            prob_val = probability_values(p);  % A appears in conditions 1-7
        end
        
        all_pref_resp = [all_pref_resp; pref_resp];
        all_pref_prob = [all_pref_prob; prob_val];
    end
    
    % Process non-preferred tone responses
    for p = 1:7
        if preferred_tone_direction(i) == 1
            % Direction 1: non-preferred tone is A
            nonpref_resp = mean(allAmat(neuron, testing_indices, p));
            prob_val = probability_values(p);  % A appears in conditions 1-7
        else
            % Direction -1: non-preferred tone is B
            nonpref_resp = mean(allBmat(neuron, testing_indices, p));
            prob_val = probability_values(p+1);  % B appears in conditions 2-8
        end
        
        all_nonpref_resp = [all_nonpref_resp; nonpref_resp];
        all_nonpref_prob = [all_nonpref_prob; prob_val];
    end
end

% Calculate correlations
[rho_om, p_om] = corr(all_om_prob, all_om_resp, 'Type', 'Spearman');
[rho_pref, p_pref] = corr(all_pref_prob, all_pref_resp, 'Type', 'Spearman');
[rho_nonpref, p_nonpref] = corr(all_nonpref_prob, all_nonpref_resp, 'Type', 'Spearman');

% Display results
fprintf('\nSpearman correlations using neuronal data:\n');
fprintf('Omission responses: r = %.4f, p = %.2e, n = %d data points\n', rho_om, p_om, length(all_om_resp));
fprintf('Preferred tone (O_P) responses: r = %.4f, p = %.2e, n = %d data points\n', rho_pref, p_pref, length(all_pref_resp));
fprintf('Non-preferred tone (O_NP) responses: r = %.4f, p = %.2e, n = %d data points\n', rho_nonpref, p_nonpref, length(all_nonpref_resp));