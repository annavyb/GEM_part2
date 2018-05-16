function [ position channame ] = readFileChanloc( filename, nelectrodes )
% reads the chanloc info and puts it into the the variables "position" and
% "channame"
%
%--------------------------------------------------------------------------
%
%  Input:
%       'filename' -[string] - [the name of the file containing the channel
%                               locations]
%  Output:
%        'position' - [ double (nelectrodes x 3)] - [ the 3d coordinates 
%                                                     of the electrodes]
%        'channame' - [cell of strings (1xnelectrodes)] - [ the channels' 
%                                                           labels]
% -------------------------------------------------------------------------
%
%  Version: V1 2017 04 29 
%  Author: VYBORNOVA Anna 
% 
% ------------------------- BEGIN CODE ------------------------------------


fid = fopen(filename);

for len = 1:nelectrodes
    position(len,:) = fscanf(fid,'%f',[3])';
    channame{len} = fscanf(fid,'%s',[1]);
end

fclose(fid);
channame{1} = '1';



end

