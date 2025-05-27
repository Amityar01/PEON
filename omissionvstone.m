%% omissionvstone.m
% This script performs the analysis and generates the pie charts corresponding
% to Figure 4C of the manuscript (Yaron et al., 2025).
% It classifies identified PEONs based on their selectivity to omissions and
% to tone presentations (Omission Preferred - OP, vs. Omission Non-Preferred - ONP).
% The script uses a split-half approach: PEONs are identified on a training
% dataset (odd trials of allommat), and their response characteristics for
% classification are evaluated on a testing dataset (even trials of allommat,
% allAmat, allBmat).
%
% Workflow:
% 1. Defines parameters: alpha (significance level), testing_indices (e.g., even trials),
%    num_probabilities. 
% 2. Prepares training_data (allommat from training_indices) and testing_data (allommat
%    from testing_indices).
% 3. Identifies PEON candidates ('PEONs_training') and their overall preferred tone
%    directions ('preferred_tone_direction_all') based on Spearman correlation
%    on the 'training_data'.
% 4. Aligns omission responses for both training ('PEON_resp_tr_al') and testing
%    ('PEON_resp_test_al') halves based on 'preferred_tone_direction_all'.
% 5. Refines the 'PEONs_training' list using a Wilcoxon signed-rank test on the
%    aligned TRAINING data ('PEON_resp_tr_al') for high probability omission bins.
%    Updates 'preferred_tone_direction' for this refined list.
% 6. For the refined PEONs, calculates their mean responses for OP tones, ONP tones,
%    and specific omission conditions (om95a0b5 - high prob OP omission,
%    om95b0a5 - high prob ONP omission, after alignment) using the TESTING_INDICES
%    applied to 'allAmat', 'allBmat', and 'allommat'.
% 7. Calls the helper function 'calculateSignificanceMasks' to categorize each
%    PEON based on its response significance (using signrank tests on trial-by-trial
%    data from testing_indices) to:
%    a) OP vs. ONP tones (categories: non-selective, OP-only, ONP-only, both).
%    b) The two types of high-probability omission conditions (categories:
%       non-selective, selective for one type, selective for the other, both).
% 8. Generates a figure with two pie charts:
%    a) Left pie: Omission selectivity (based on sig_masks_omission).
%    b) Right pie: For omission-selective PEONs, their tone response selectivity
%       (based on sig_masks_tone).
% 9. Saves the figure as 'panel_C.pdf'.
%
% Inputs (expected to be in the MATLAB workspace):
%   - allommat: 3D matrix (neurons x 50 trials x 8 probability_conditions) of
%               trial-by-trial omission responses (from 'allmats.m').
%   - allAmat:  3D matrix (neurons x 50 trials x 7 tone_A_conditions) (from 'allmats.m').
%   - allBmat:  3D matrix (neurons x 50 trials x 7 tone_B_conditions) (from 'allmats.m').
%  
%                  MUST be defined in the workspace before running this script.
%
% Outputs:
%   - 'panel_C.pdf': Saved figure with the two pie charts.
%   - Console output: "Panel C saved as panel_C.pdf".
%   - Workspace variables created/modified:
%     - PEONs_training: Indices of PEONs as defined by this script's criteria for Fig 4C.
%     - preferred_tone_direction: Preferred direction for the PEONs in this script's list.
%     - sig_masks_tone, sig_masks_omission: Structures with logical masks for
%       neuron classifications used in the pie charts.
%     - And other intermediate variables (rho_train, p_train, Onp, Op, etc.).
%
% Dependencies:
%   - Standard MATLAB plotting and statistics functions (e.g., 'corr', 'signrank', 'pie').
%   - Helper function 'calculateSignificanceMasks' (defined within this script).
%
% Author: Amit Yaron

alpha = 0.05;                        % significance level
testing_indices = 2:2:50;            % EVEN trials for test set
num_probabilities = 8;

% -------------------------------------------------------------
% 1)  |rho| & p-values on TRAINING vs TEST  halves
% -------------------------------------------------------------

training_data = allommat(:, training_indices, :);
testing_data  = allommat(:, testing_indices, :);

num_neurons  = size(allommat,1);

prob_labels_train = repmat(1:num_probabilities , 1, numel(training_indices))';
prob_labels_test  = repmat(1:num_probabilities , 1, numel(testing_indices))';

training_flat = reshape(permute(training_data,[1 3 2]), num_neurons, [])';
testing_flat  = reshape(permute(testing_data, [1 3 2]), num_neurons, [])';

rho_train = zeros(num_neurons,1);
p_train   = zeros(num_neurons,1);
rho_test  = zeros(num_neurons,1);
for n = 1:num_neurons
    [rho_train(n), p_train(n)] = corr(prob_labels_train, training_flat(:,n), 'type','Spearman');
    [rho_test(n) , ~        ] = corr(prob_labels_train, testing_flat(:,n) , 'type','Spearman');
end

PEONs_training = find(p_train < alpha);               % PEON candidates
preferred_tone_direction_all = sign(rho_train);       % +1 or –1 for every neuron
preferred_tone_direction      = preferred_tone_direction_all(PEONs_training);

% -------------------------------------------------------------
% 2) Align omission responses (flip if rho<0)
% -------------------------------------------------------------
PEON_resp_test = testing_data(PEONs_training,:,:);    % [nPEON x 25 x 8]
PEON_resp_tr   = training_data;                       % all neurons

% flip in TEST half
PEON_resp_test_al = zeros(size(PEON_resp_test));
for i = 1:numel(PEONs_training)
    tmp = squeeze(PEON_resp_test(i,:,:));             % 25 x 8
    if preferred_tone_direction(i)==-1, tmp = fliplr(tmp); end
    PEON_resp_test_al(i,:,:) = tmp;
end

% flip in TRAIN half (for Wilcoxon)
PEON_resp_tr_al = zeros(size(PEON_resp_tr));
for i = 1:num_neurons
    tmp = squeeze(PEON_resp_tr(i,:,:));               % 25 x 8
    if preferred_tone_direction_all(i)==-1, tmp = fliplr(tmp); end
    PEON_resp_tr_al(i,:,:) = tmp;
end

% -------------------------------------------------------------
% 3) Wilcoxon on high-prob bins (5:8)  → refined PEONs
% -------------------------------------------------------------
high_prob_bins = 5:8;
is_resp_sig = false(num_neurons,1);
for i = 1:num_neurons
    block = PEON_resp_tr_al(i,:,high_prob_bins);
    [pWil,~] = signrank(block(:),0);
    is_resp_sig(i) = pWil < alpha;
end
valid_mask = is_resp_sig(PEONs_training);             % logical for PEONs_training
PEONs_training = PEONs_training(valid_mask);
preferred_tone_direction = preferred_tone_direction(valid_mask);
PEON_resp_test_al = PEON_resp_test_al(valid_mask,:,:);

% -------------------------------------------------------------
% 4)  Compute Op/Onp & omission 95/5 responses (EVEN trials)
% -------------------------------------------------------------
Onp = zeros(numel(PEONs_training),1);
Op  = zeros(numel(PEONs_training),1);
om95a0b5 = zeros(numel(PEONs_training),1);
om95b0a5 = zeros(numel(PEONs_training),1);

for i = 1:numel(PEONs_training)
    idx = PEONs_training(i);
    if preferred_tone_direction(i)==1
        Onp(i)       = mean(allAmat(idx, testing_indices, 7), 'all');
        Op(i)        = mean(allBmat(idx, testing_indices, 1), 'all');
        om95b0a5(i)  = mean(allommat(idx, testing_indices, 1), 'all');
        om95a0b5(i)  = mean(allommat(idx, testing_indices, 8), 'all');
    else
        Op(i)        = mean(allAmat(idx, testing_indices, 7), 'all');
        Onp(i)       = mean(allBmat(idx, testing_indices, 1), 'all');
        om95a0b5(i)  = mean(allommat(idx, testing_indices, 1), 'all');
        om95b0a5(i)  = mean(allommat(idx, testing_indices, 8), 'all');
    end
end

% -------------------------------------------------------------
% 5) Significance masks helper
% -------------------------------------------------------------
[sig_masks_tone, sig_masks_omission] = calculateSignificanceMasks( ...
        Onp, Op, om95a0b5, om95b0a5, ...
        allAmat, allBmat, allommat, ...
        PEONs_training, preferred_tone_direction, testing_indices);

% -------------------------------------------------------------
% 6) ===  FIGURE C  ===  (pie charts exactly as screenshot)
% -------------------------------------------------------------
Cfig = figure('Position',[100 100 1000 400]);

% -- LEFT pie  (59 / 31)
subplot(1,2,1)
om_sel  = sig_masks_omission.xonly;
om_non  = sig_masks_omission.both;
valsL   = [sum(om_non)  sum(om_sel)];          % [31 59]
labsL   = {sprintf('Omission\nnon-selective\n(%d)',valsL(1)), ...
           sprintf('Omission\nselective\n(%d)',valsL(2))};
colsL   = [0.55 0.8 0.3; 0.2 0.6 1];

hL = pie(valsL);  colormap(gca,colsL);  delete(hL(2:2:end));
text([-0.5 0.5],[0.5 -0.5], labsL,'Horiz','center','FontSize',12);
title({'\bfC','PEON_{odd}','in EVEN trials'},'FontSize',14);

% -- RIGHT pie (3 / 6 / 41 / 9)
subplot(1,2,2)
xonly = sig_masks_omission.xonly;             % 59 neurons
valsR = [ ...
    sum(sig_masks_tone.xonly & xonly),  ...    % Op-only  (3)
    sum(sig_masks_tone.non   & xonly),  ...    % None     (6)
    sum(sig_masks_tone.both  & xonly),  ...    % Both     (41)
    sum(sig_masks_tone.yonly & xonly)];        % Onp-only (9)

labsR = {sprintf('O_P-only\n(%d)',valsR(1)), ...
         sprintf('None\n(%d)',valsR(2)), ...
         sprintf('Both\n(%d)',valsR(3)), ...
         sprintf('O_{NP}-only\n(%d)',valsR(4))};

colsR = [1 0.6 0.2; 1 1 0.5; 0.6 0 0; 1 0.4 0.1];

hR = pie(valsR);  colormap(gca,colsR);  delete(hR(2:2:end));
text([-0.7 0 0.5 0],[0.5 0.8 0 -0.8],labsR,'Horiz','center','FontSize',12);
title('Omission selective neurons','FontSize',12);

% -- Arrow between pies
annotation('arrow',[0.47 0.53],[0.5 0.5],'LineWidth',2);

set(Cfig,'PaperPositionMode','auto');
print(Cfig,'panel_C.pdf','-dpdf','-painters','-r600');

disp('Panel C saved as panel_C.pdf');

%% ------------------------------------------------------------------------
% Helper: Calculate significance masks  (UNCHANGED LOGIC)
% ------------------------------------------------------------------------
function [sig_masks_tone, sig_masks_omission] = calculateSignificanceMasks( ...
    Onp, Op, om95a0b5, om95b0a5, ...
    allAmat, allBmat, allommat, ...
    PEONs_training, prefDir, idx)
% tone
sig_Onp = false(size(Onp)); sig_Op  = false(size(Op));
for k = 1:numel(Onp)
    if prefDir(k)==1
        reps_Onp = squeeze(allAmat(PEONs_training(k),idx,1:7));
        reps_Op  = squeeze(allBmat(PEONs_training(k),idx,1:7));
    else
        reps_Onp = squeeze(allBmat(PEONs_training(k),idx,1:7));
        reps_Op  = squeeze(allAmat(PEONs_training(k),idx,1:7));
    end
    sig_Onp(k) = signrank(reps_Onp(:),0,'tail','right') < 0.05 & mean(reps_Onp,'all')>0;
    sig_Op(k)  = signrank(reps_Op(:) ,0,'tail','right') < 0.05 & mean(reps_Op ,'all')>0;
end
sig_masks_tone.non   = ~sig_Onp & ~sig_Op;
sig_masks_tone.xonly =  sig_Op & ~sig_Onp;
sig_masks_tone.yonly =  sig_Onp & ~sig_Op;
sig_masks_tone.both  =  sig_Op &  sig_Onp;

% omission
sig_a = false(size(om95a0b5)); sig_b = false(size(om95b0a5));
for k = 1:numel(om95a0b5)
    if prefDir(k)==-1
        reps_a = squeeze(allommat(PEONs_training(k),idx,1:4));
        reps_b = squeeze(allommat(PEONs_training(k),idx,5:8));
    else
        reps_a = squeeze(allommat(PEONs_training(k),idx,5:8));
        reps_b = squeeze(allommat(PEONs_training(k),idx,1:4));
    end
    sig_a(k) = signrank(reps_a(:),0,'tail','right')<0.05 & mean(reps_a,'all')>0;
    sig_b(k) = signrank(reps_b(:),0,'tail','right')<0.05 & mean(reps_b,'all')>0;
end
sig_masks_omission.non   = ~sig_a & ~sig_b;
sig_masks_omission.xonly =  sig_a & ~sig_b;
sig_masks_omission.yonly = ~sig_a &  sig_b;
sig_masks_omission.both  =  sig_a &  sig_b;
end
