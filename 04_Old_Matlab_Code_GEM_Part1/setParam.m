function [Params] = setParam(ShowParam)
%This function sets the principle to be used in function
%SeizureDataDiscoveringV_04_29
if nargin == 0
    ShowParam = false; 
end
Params.FeatureRanking.Method = 'RANK_SUM';% Can take values :
% 1) 'RANK_SUM';
% 2) 'FISHER';
% 3) 'TTEST';
Params.FeatureRanking.Values = 'MEAN'; % Can take values :
% 1) 'ALL_CHANNELS';
% 2) 'MEAN'
% 3) 'VARIANCE'
% 4) 'MEAN_AND_VARIANCE'
Params.BadChannels = 'INTERPOLATE'; % Can take the values :
% 1) 'REMOVE'
% 2) 'INTERPOLATE'

% --- FILTERING ----

Params.Filtering.FiltOrder = 8; % filter order 
Params.Filtering.LowCf = 1; % lower cut-off frequency
Params.Filtering.HighCf = 30; % higher cut-off frequency 

% --- DURATIONS OF THE SEGMENTS ( in seconds ) ---

Params.Duration.Interictal = 5*60; % non-seizure data window
Params.Duration.Ictal = 2*60; % seizure data window (goes after the seizure 
%                               marker)
Params.Duration.PreIctal = 3*60; % pre-seizure data window (goes just before 
%                                  marker)

%--- GFT PARAMETERS ---

Params.Gft.R = 80; % the Radius in mm defining the neighborhood for each 
%                    node of the graph 

Params.PreliminaryStatistics = true; 
Params.PlotAll = false; 

if ShowParam
    disp('PARAMETERS:');
    disp(['T-test method:',Params.FeatureRanking.Method]);
    disp(['T-test computed on: ', Params.FeatureRanking.Values]);
    disp(['Bad Channels handling: ', Params.BadChannels ]);
    disp(['Filtering between [', num2str(Params.Filtering.LowCf), ',',...
        num2str(Params.Filtering.HighCf), ']', ' order ', ...
        num2str(Params.Filtering.FiltOrder)]); 
    disp(['Duration window before seizure: ', ...
        num2str(Params.Duration.Interictal), ' [s]']); 
    disp(['Duration window pre-seizure: ', ...
        num2str(Params.Duration.PreIctal), ' [s]']);
    disp(['Duration window during seizure: ', ...
        num2str(Params.Duration.Ictal), ' [s]']); 
    disp(['R - radius of the neighborhood: ', ...
        num2str(Params.Gft.R), ' [mm]']); 
    
end

end

