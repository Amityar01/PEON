
%% batchall.m
% Master script for batch processing of neuronal data from multiple datasets.
% Iterates through datasets defined in 'tables.mat', loads raw data and
% corresponding baseline data. For each dataset, it calls
% 'maketests1214NosideAllAud.m' to perform statistical tests and calculate
% mean responses/latencies.
%
% Workflow:
% 1. Loads 'tables.mat' (dataset identifiers and metadata).
% 2. For each dataset:
%    a. Loads 'res[dataset_id].mat' (raw responses AND trial condition indices like inddss, inddds).
%    b. Loads 'bl[dataset_id].mat' (baseline data).
%    c. Calls 'maketests1214NosideAllAud.m' with 'res' and 'bl'. The 'ttt' variable is
%       incremented here to keep a running count of processed neurons across datasets.
%    d. Stores outputs from 'maketests1214NosideAllAud.m' (e.g., meanom, test)
%       into dynamically named workspace variables (e.g., meanomC1, testC1).
%
% Inputs (expected in MATLAB path or workspace):
%   - tables.mat: Cell array defining dataset identifiers and metadata.
%                 Example row: {'C1', ToneAfreq, ToneBfreq, RatID, PenetrationID, LocationStr}.
%   - Raw data files for each dataset ID in 'tables.mat':
%     - 'res[dataset_id].mat': Should contain 'res' (trials x timepoints x neurons data)
%                              and trial index variables (e.g., indd5, indd55A, etc.,
%                              which are then used by 'maketests1214NosideAllAud.m'
%                              to define 'inddss', 'inddds').
%     - 'bl[dataset_id].mat': Should contain 'bl' (baseline data, similar structure to 'res').
%
% Outputs:
%   - Dynamically named variables in the MATLAB workspace for each dataset,
%     e.g., meanomC1, testC1, chC1 (channel/depth info). These are then typically
%     used by 'createtable.m' and subsequently 'allmats.m'.
%   - Cell arrays 'allom', 'allst', 'alldv', 'allind' are accumulated across datasets
%     in the workspace, containing trial-by-trial responses and indices.
%
% Dependencies:
%   - maketests1214NosideAllAud.m
%
% Author: Amit Yaron
clear all
close all
sur=10*[260,220,305,222,284,246,300,270,170,210,260,280,280,280,340,220,240,384+40]-520;
RR=[];
LL=[];
RRd=[]; 
LLd=[];
ttt=0;
load tables
for tt =1:18
    eval(['load(''res',tables{1,tt },'.mat'')']);
    eval(['load(''bl',tables{1,tt },'.mat'')']);
   % eval(['bl=bl',tables{1,tt },';'])
    maketests1214NosideAllAud       ;
    eval(['meanom',tables{1,tt },'=meanom;']);
     eval(['ch',tables{1,tt },'=((sur(tt)-ch));']);
    eval(['meanst',tables{1,tt },'=meanst;']);
    eval(['meandv',tables{1,tt },'=meandv;']);
    eval(['latom',tables{1,tt },'=latom;']);
    eval(['latst',tables{1,tt },'=latst;']);
    eval(['latdv',tables{1,tt },'=latdv;']);
    eval(['test',tables{1,tt },'=test;']);
    clear rr
    clear ll
    for i=1:length(nu)
        if side(i)==1
            rr(:,i)=squeeze(nanmean((res(indd5(1:end)+1000,1:600,nu(i))),1));
            ll(:,i)=squeeze(nanmean((res(indd5(1:end),1:600,nu(i))),1));
            rrd(:,i)=squeeze(nanmean((res(indd55A-4000,1:600,nu(i))),1));
            lld(:,i)=squeeze(nanmean((res(indd55B-4000,1:600,nu(i))),1));
        else
            rr(:,i)=squeeze(nanmean((res(indd5(1:end),1:600,nu(i))),1));
            ll(:,i)=squeeze(nanmean((res(indd5(1:end)+1000,1:600,nu(i))),1));
            lld(:,i)=squeeze(nanmean((res(indd55A-4000,1:600,nu(i))),1));
            rrd(:,i)=squeeze(nanmean((res(indd55B-4000,1:600,nu(i))),1));
        end
    end
    RR=[RR, rr];
    LL =[LL, ll];
    RRd=[RRd, rrd];
    LLd =[LLd, lld];
end
createtable 