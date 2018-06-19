function [ psorted, ind, indSignif ] = calculateStatistics(Metric,Labels, alpha, Label1, Label2, method)

if size(Metric,2) ~= 204
    warning(' number of features different from 204')
end
Ind1 = find( Labels == Label1); % underrepresented class
Ind2 = find(Labels == Label2); % overrepresented class

Metric1 = Metric(Ind1, :);
Metric2 = Metric(Ind2, :);
l = randperm(length(Ind2),length(Ind1));

if nargin == 6
    for iF = 1:size(Metric,2)
        if strcmp(method, 'ranksum')
            [p(iF), h(iF)] = ranksum(Metric1(:,iF), Metric2(:,iF));
        elseif strcmp(method, 'ttest2')
            [p(iF), h(iF)] = ttest2(Metric1(:,iF), Metric2(:,iF));
        elseif strcmp(method, 'signrank')
            [p(iF), h(iF)] = signrank(Metric1(:,iF), Metric2(l,iF));
        elseif strcmp(method, 'ttest')
            [p(iF), h(iF)] = signrank(Metric1(:,iF), Metric2(l,iF));
        end
    end
else
    
    for iF = 1:size(Metric,2)
        
        [p(1,iF), ~] = ranksum(Metric1(:,iF), Metric2(:,iF));
        [~, p(2,iF)] = ttest2(Metric1(:,iF), Metric2(:,iF));
        [p(3,iF), ~] = signrank(Metric1(:,iF), Metric2(l,iF));
        [~, p(4,iF)] = ttest(Metric1(:,iF), Metric2(l,iF));
        
    end
end

    [psorted, ind] = sort(p,2); 
    
    % find significant indeces alpha = 0.01
    
    for iP = 1: size(psorted,1)
        indSignif{iP} = find((p(iP,:) < alpha) == 1); 
    end

end

