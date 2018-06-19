function [ err ] = classerr( l1, l2 )
ind = 0; 
for il = 1:length(l1)
    if l1(il) ~= l2(il)
       ind = ind+1; 
    end
end
err = 100*ind/length(l1); 
end

