%% maketests1214NosideAllAud.m
% Processes neuronal data for a SINGLE dataset/penetration.
% Called by 'batchall.m'. This script performs statistical tests and extracts
% response metrics for omission, standard, and deviant tone conditions.
%
% Core Functions:
% 1. Identifies active neurons ('nu') based on activity in 'res'.
% 2. Uses trial index variables (e.g., indd5, indd55A, loaded from the
%    res[dataset_id].mat file by batchall.m and thus available in the
%    workspace) to define 'inddss' and 'inddds' for different conditions.
% 3. For active neurons and each condition:
%    a. Performs Wilcoxon signed-rank tests (response window vs. baseline window).
%    b. Calculates mean baseline-subtracted responses: 'meanom', 'meanst', 'meandv'.
%    c. Extracts trial-by-trial responses into 'test.allom', 'test.allst', 'test.alldv'.
%    d. Estimates response latencies: 'latom', 'latst', 'latdv'.
% 4. Stores results in output variables.
%
% Inputs (passed from or available in the workspace of 'batchall.m'):
%   - res: 3D matrix (trials x timepoints x neurons) of neuronal responses for the current dataset.
%          This .mat file (loaded by batchall.m) must also contain the trial index
%          variables (e.g., indd5, indd55A) used to create inddss/inddds.
%   - bl: 3D matrix (trials x timepoints x neurons) of baseline activity.
%   - ttt: (Optional, if used for context from batchall.m) Scalar, running total of neurons processed
%          before this dataset, used for indexing into global 'allom', 'allst', etc.
%          This script itself primarily processes the current dataset's neurons.
%
% Outputs (variables created/modified in the caller's workspace, i.e., batchall.m):
%   - meanom: Matrix (active_neurons x conditions) of mean omission responses.
%   - meanst: Matrix (active_neurons x conditions) of mean standard tone responses.
%   - meandv: Matrix (active_neurons x conditions) of mean deviant tone responses.
%   - latom, latst, latdv: Matrices for response latencies.
%   - test: Structure containing detailed results, including 'test.allom',
%           'test.allst', 'test.alldv' (cell arrays of trial-by-trial responses),
%           'test.TT' (p-values from signrank for various response types).
%   - T: Matrix (conditions x active_neurons) containing primary signrank p-values for omissions.
%   - nu: Vector of indices for active neurons within the current dataset.
%   - (Modifies global-like cell arrays: allom, allst, alldv, allind by appending
%     trial-by-trial data using ttt+t as row index).
%
% Note: Specific time windows for baseline and response are hardcoded.
%
% Author: Amit Yaron

%%
clear test
clear T
clear TT*
nu=find(((squeeze(sum(sum(res(:,:,:)))))>200)');
inddss={indd5+1000,indd55+2000,indd105+4000,indd205+6000,indd205+7000,indd105+5000,indd55+3000,indd5};
inddds={indd55A-4000,indd55A-4000,indd105A-4000,indd205A-4000,indd205B-4000,indd105B-4000,indd55B-4000,indd55B-4000};
T=zeros(size(nu));
test.totalneurons=size(res,3);
test.activeneurons=length(nu);

clear A
clear omnu
clear Anu
clear Bnu

t=0;
for ii=nu;
    t=t+1;
    for jj=[1:8]
        [TT(round(jj),t),~,~]=(signrank(squeeze(mean(res(inddss{jj},306:420,ii),2)),squeeze(mean(res(inddss{jj}(:),276:305,ii),2)),'tail','right','alpha',0.05)); % 1-8 omission
        [TT(round(jj)+8,t),~,~]=(signrank(squeeze(mean(res(inddss{jj},(306:420)-150,ii),2)),squeeze(mean(res(inddss{jj}(:),(276:305)-150,ii),2)),'tail','right','alpha',0.05)); % 9-16 standard
    end
    for jj=[2:7]
        [TT(jj+15,t),~,~]=(signrank(squeeze(mean(res(inddds{jj},306:420,ii),2)),squeeze(mean(res(inddds{jj}(:),276:305,ii),2)),'tail','right','alpha',0.05));% 17-22 deviant
      
    end
    
end
alpha=0.05;
%auditory neurons
num_conditions=22;
for ii=1:(size(TT,2))
     A(ii)=1;
%     p_values = TT(:,ii)
%     sorted_p_values = sort(p_values);
%     critical_values = (1:num_conditions) * (alpha / num_conditions);
%     % Find the largest p-value that is smaller than its critical value
%     max_p_value_index = find(sorted_p_values' <= critical_values, 1, 'last');
%     % Check if any condition is significant
%     A(ii) = false;
%     if ~isempty(max_p_value_index)
%         A(ii) = any(p_values <= sorted_p_values(max_p_value_index));
%     else
%         A(ii)=0;
%     end
end

% omission neurons
num_conditions=8;
for ii=1:(size(TT,2))
    p_values = TT(1:8,ii);
    sorted_p_values = sort(p_values);
    critical_values = (1:num_conditions) * (alpha / num_conditions);
    % Find the largest p-value that is smaller than its critical value
    max_p_value_index = find(sorted_p_values' <= critical_values, 1, 'last');
    % Check if any condition is significant
    Onu(ii) = false;
    if ~isempty(max_p_value_index)
        Onu(ii) = any(p_values <= sorted_p_values(max_p_value_index));
    else
        Onu(ii) = 0;
    end
end

% A neurons

%Anu=TT(22,:)<0.05;
num_conditions=7;
for ii=1:(size(TT,2))                      
    p_values = TT([9 10 11 12 20 21  22],ii);
    sorted_p_values = sort(p_values);
    critical_values = (1:num_conditions) * (alpha / num_conditions);
    % Find the largest p-value that is smaller than its critical value
    max_p_value_index = find(sorted_p_values' <= critical_values, 1, 'last');
    % Check if any condition is significant
    Anu(ii) = false;
    if ~isempty(max_p_value_index)
        Anu(ii) = any(p_values <= sorted_p_values(max_p_value_index));
    else
        Anu(ii) = 0;
    end
end

% B neurons
Bnu=TT(17,:)<0.05;
num_conditions=7;
for ii=1:(size(TT,2))   
    p_values = TT([17 18 19 13 14 15 16 ],ii);
    sorted_p_values = sort(p_values);
    critical_values = (1:num_conditions) * (alpha / num_conditions);
    % Find the largest p-value that is smaller than its critical value
    max_p_value_index = find(sorted_p_values' <= critical_values, 1, 'last');
    % Check if any condition is significant
    Bnu(ii) = false;
    if ~isempty(max_p_value_index)
        Bnu(ii) = any(p_values <= sorted_p_values(max_p_value_index));
    else
        Bnu(ii) = 0;
    end
end

nu=nu(A>0);
test.auditoryneurons=sum(A>0);
%test.omissiont=omissiont;
test.omissionall=sum(TT(1,(A>0))<0.025|TT(2,(A>0))<0.025);
test.allp=test.omissionall/test.auditoryneurons;
omboth=(TT(1,(A>0))<0.05)&( TT(8,(A>0))<0.05);
omnu=Onu;
%nu=nu(omnu>0);

TT=TT(:,A>0);
Test.TT=TT;
%%
%filename = [Date,'continuous\Neuropix-PXI-105.0\output\cluster_info.tsv'];
%T = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');
for i=1:length(nu); ch(i)=cluster_depths(cluster_ids==uni(nu(i))) ;end
clear meanom
clear meanst
clear meandv
clear latom
clear latst
clear latdv
% inddss={indd5+1000,indd55+6000,indd105+8000,indd205+10000,indd205+11000,indd105+9000,indd55+7000,indd5};
% indsss={inds5A,inds55A,inds105A,inds205A,inds205B,inds105B,inds55B,inds5B};
% inddds={1,indd55A,indd105A,indd205A,indd205B,indd105B,indd55B,1};
inddss={indd5+1000,indd55+2000,indd105+4000,indo205A-4000,indo205B-4000,indd105+5000,indd55+3000,indd5};
inddds={indd55A-4000,indd55A-4000,indd105A-4000,indd205A-4000,indd205B-4000,indd105B-4000,indd55B-4000,indd55B-4000};
indsss={inds5B,inds55B-4000,inds105B-4000,inds205B'-4000,inds205A'-4000,inds105A-4000,inds55A-4000,inds5A};
for i=1:8
indsss{i}=setdiff(indsss{i},inddss{i}+1);
end

t=0;
clear r
clear l
for ii=nu
    t=t+1;
    l(t)=1;%(mean(smoothdata(mean(res([inddss{1}],306:420,ii),2),'gaussian',1))-mean(mean(res([inddss{1}],276:305,ii))));
    r(t)=0;%(mean(smoothdata(mean(res([inddss{8}],306:420,ii),2),'gaussian',1))-mean(mean(res([inddss{8}],276:305,ii))));

end
side=(r>l);
for jj=1:8
    t=0;
    for  ii=nu
        t=t+1;
        
        if side(t)==0
             % Omission response
            baseline_mean = mean(mean(res([inddss{jj}], 276:305, ii), 2));
            baseline_peak = max(smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20));%std(mean(res([inddss{jj}], 276:305, ii), 2));
            threshold = 0*baseline_mean + baseline_peak*0.5; % Define the threshold
            response_activity = smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20);
            latency_indices = find(response_activity > threshold, 1);
            if ~isempty(latency_indices)
                latom(t, jj) = latency_indices+5; % Adjust index to match time window
            else
                latom(t, jj) = NaN;
            end
            meanom(t,jj)=(mean(smoothdata(mean(res([inddss{jj}],306:420,ii),2),'gaussian',1))-mean(mean(res([inddss{jj}],276:305,ii))));
            allom{ttt+t,jj}=mean(res([inddss{jj}],306:420,ii),2)'-mean(res([inddss{jj}],276:305,ii),2)';
            allind{ttt+t,jj}=inddss{jj};
             % Standard response
            baseline_mean = mean(mean(res([inddss{jj}], [276:305] - 150, ii)));  
            baseline_peak = max(smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20));%std(nanmean(res([inddss{jj}], [276:305] - 150, ii), 2));
            threshold = 0*baseline_mean +  baseline_peak*0.5;

            response_activity = smoothdata(mean(res([inddss{jj}], (306:420) - 150, ii), 2),'gaussian',20); 
            latency_indices = find(response_activity > threshold, 1);
            if ~isempty(latency_indices)
                latst(t, jj) = latency_indices+5; % Adjust index to match time window
            else
                latst(t, jj) = NaN;
            end
            meanst(t,jj)=(mean(smoothdata(mean(res([inddss{jj}],[306:420]-150,ii),2),'gaussian',1))-1*mean(mean(res([inddss{jj}],[276:305]-150,ii))));
            allst{ttt+t,jj}=mean(res([indsss{jj}],[306:420],ii),2)'-1*mean(res([indsss{jj}],[276:305],ii),2)';
            
            % Deviant response
            if jj <= 7 % Only calculate for valid deviant indices
                baseline_mean = mean(mean(res([inddds{jj}], 276:305, ii), 2));
                baseline_peak = max(smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20));%std(mean(res([inddds{jj}], 276:305, ii), 2));
                threshold = 0*baseline_mean + baseline_peak*0.5;

                response_activity = smoothdata(mean(res([inddds{jj}], 306:420, ii), 2),'gaussian',20);
                latency_indices = find(response_activity > threshold, 1);
                if ~isempty(latency_indices)
                    latdv(t, jj) = latency_indices+5; % Adjust index to match time window
                else
                    latdv(t, jj) = NaN;
                end
            end
            meandv(t,jj)=(mean(smoothdata(mean(res([inddds{jj}],306:420,ii),2),'gaussian',1))-1*mean(mean(res([inddds{jj}],276:305,ii))));
            alldv{ttt+t,jj}=mean(res([inddds{jj}],[306:420],ii),2)'-1*mean(res([inddds{jj}],[276:305],ii),2)';
           
           
        else
              % Omission response
            baseline_mean = mean(mean(res([inddss{jj}], 276:305, ii), 2));
            baseline_peak = max(smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20));%std(mean(res([inddss{jj}], 276:305, ii), 2));
            threshold = 0*baseline_mean +  baseline_peak*0.5;

            response_activity = smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20);
            latency_indices = find(response_activity > threshold, 1);
            if ~isempty(latency_indices)
                latom(t, 9 - jj) = latency_indices+5;
            else
                latom(t, 9 - jj) = NaN;
            end

            % Standard response
            baseline_mean = mean(mean(res([inddss{jj}], (276:305) - 150, ii), 2));
            baseline_peak = max(smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20));%std(mean(res([inddss{jj}], (276:305) - 150, ii), 2));
            threshold = 0*baseline_mean +  baseline_peak*0.5;

            response_activity = smoothdata(mean(res([inddss{jj}], (306:420) - 150, ii), 2),'gaussian',20);
            latency_indices = find(response_activity > threshold, 1);
            if ~isempty(latency_indices)
                latst(t, 9 - jj) = latency_indices+5;
            else
                latst(t, 9 - jj) = NaN;
            end

            % Deviant response
            if jj <= 7 % Only calculate for valid deviant indices
                baseline_mean = mean(mean(res([inddds{jj}], 276:305, ii), 2));
                baseline_peak = max(smoothdata(mean(res([inddss{jj}], 306:420, ii), 2),'gaussian',20));%std(mean(res([inddds{jj}], 276:305, ii), 2));
                threshold = 0*baseline_mean +  baseline_peak*0.5;

                response_activity = smoothdata(mean(res([inddds{jj}], 306:420, ii), 2),'gaussian',20);
                latency_indices = find(response_activity > threshold, 1);
                if ~isempty(latency_indices)
                    latdv(t, 9 - jj) = latency_indices+5;
                else
                    latdv(t, 9 - jj) = NaN;
                end
            end
            meanom(t,9-jj)=(mean(smoothdata(mean(res([inddss{jj}],306:420,ii),2),'gaussian',1))-mean(mean(res([inddss{jj}],276:305,ii))));
            allom{ttt+t,9-jj}=mean(res([inddss{jj}],306:420,ii),2)'-mean(res([inddss{jj}],276:305,ii),2)';
            allind{ttt+t,9-jj}=inddss{jj};
            meanst(t,9-jj)=(mean(smoothdata(mean(res([inddss{jj}],[306:420]-150,ii),2),'gaussian',1))-1*mean(mean(res([inddss{jj}],[276:305]-150,ii))));
            allst{ttt+t,9-jj}=mean(res([indsss{jj}],[306:420],ii),2)'-1*mean(res([indsss{jj}],[276:305],ii),2)';
            meandv(t,9-jj)=(mean(smoothdata(mean(res([inddds{jj}],306:420,ii),2),'gaussian',1))-1*mean(mean(res([inddds{jj}],276:305,ii))));
            alldv{ttt+t,9-jj}=mean(res([inddds{jj}],[306:420],ii),2)'-1*mean(res([inddds{jj}],[276:305],ii),2)';
             temp=Anu(t);
             Anu(t)= Bnu(t);
            Bnu(t)=temp;
        end
    end
end
ttt=ttt+t;
sel=find(~isinf(mean(meanom')) & ~isnan(mean(meanom')) & ~(mean(meanom')==0));
nu=nu(sel);
meanom=meanom(sel,:) ;
meandv=meandv(sel,:) ;
meanst=meanst(sel,:) ;
latom=latom(sel,:) ;
latdv=latdv(sel,:) ;
latst=latst(sel,:) ;
ch=ch(sel);
side=side(sel);
sch=sort(ch);
meandv(:,[1,8])=nan;
latdv(:,[1,8])=nan;
test.omnu=omnu(sel);
test.Anu=Anu(sel);
test.Bnu=Bnu(sel);
test.omboth=omboth(sel);
test.side=side;
test.omdepth=((sur(tt)-ch(test.omnu>0)));
test.TT=TT(:,sel);