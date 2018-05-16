function [ind_reorder] = convertEGItoChanlocs(chanlocs_labels,egi_labels)
% Finds the mapping from the EGI xyz channel order to the readlocs (eeglab 
% function) channel order 
%--------------------------------------------------------------------------
%
%  Input:
%       'chanlocs_labels' -[cell array] - [the labels of the channels in the 
%           order provided by readlocs function of the eeglab software]
%       'egi_labels' - [cell array] - [the labels of the channels in the 
%           order provided by egi xyz file] 
%  Output:
%        'index_reorder' - [double array] - the reordering of the indeces
%        to go from EGI configuration to the readlocs configuration
% -------------------------------------------------------------------------
%
%  Version: V1 2018 04 26 
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

indSame = 1;
for iChL = 1:length(chanlocs_labels)
    for iChO = 1:length(egi_labels)
        if strcmp(chanlocs_labels(iChL), egi_labels(iChO))
            ind_reorder(iChL) = iChO;
        end
    end
end

end

