fid = fopen('EGI 204.Geneva Average 13.10-10.xyz');
    header = fscanf(fid,'%f',[2,1]);
 
    for len = 1:header(1)
        position(len,:) = fscanf(fid,'%f',[3])';
        channame{len} = fscanf(fid,'%s',[1]);
    end
 
    fclose(fid);
