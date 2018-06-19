function [ err ] = class_Spec_err( l1, l2 )
% l1 - true labels


classes = unique(l1);
for i = 1:length(classes)
   totInd(i) = length(find(l1 == classes(i))) ; 
end
ind = zeros(1,2);

for iCl = 1:length(classes)
for il = 1:length(l1)
    
    if l1(il) == classes(iCl)
        
        if l1(il) ~= l2(il)
            ind(iCl) = ind(iCl)+1;
            
        end
    end
end
end

err = 100*ind./totInd;
end

