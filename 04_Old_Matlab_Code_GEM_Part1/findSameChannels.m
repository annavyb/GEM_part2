function [ IndSameChan, Labels ] = findSameChannels(chanlocs1, chanlocs2, nEl1, nEl2)
% finds and returns the indeces of the corresponding channels when going from 257
% to 204 channels 
%
%--------------------------------------------------------------------------
%
%  Input:
%       'chanloc1' -[string] - [the name of the file containing the channel
%                               locations]
%       'chanloc2' -[string] - [the name of the file containing the channel
%                               locations]
%  Output:
%        'IndSameChannels' - [ double] - [ the indeces of corresponding 
%                                             channels]
%
% -------------------------------------------------------------------------
%
%  Version: V1 2017 04 29 
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------

if (nEl1 > nEl2) 
    
  chanlocs257 = chanlocs1; 
  chanlocs204 = chanlocs2; 
  nEl257 = nEl1; 
  nEl204 = nEl2; 
  
elseif (nEl1 < nEl2)
    
  chanlocs257 = chanlocs2; 
  chanlocs204 = chanlocs1; 
  nEl257 = nEl2; 
  nEl204 = nEl1; 
  
end

    [position257, channame257] = readFileChanloc(chanlocs257, nEl257); 
    [position204, channame204] = readFileChanloc(chanlocs204, nEl204); 
    
     indSame = 1;

    
    for iName257 = 1:length(channame257)
       for iName204 = 1:length(channame204)
          if strcmp(channame257{iName257}, channame204{iName204})
              IndSameChan(indSame) = iName257;
              Labels{indSame} = channame204{iName204}; 
              indSame = indSame+1; 
          end
       end
    end
   
    
end

