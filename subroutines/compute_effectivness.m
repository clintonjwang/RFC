function [AUC,sensitivity,specificity,precision,accuracy,DSC]...
    = compute_effectivness(scores,labels)

%loop through labels
for c=1:3
    
    %check if patient has corresponding tissue class
    if(length(find(labels==c))>0)
        lab = double(labels==c);
        [~,~,~,AUC(c,1)] = perfcurve(lab,scores(:,c+1),1);
        
        if(c==1)
            classification = (scores(:,2)>scores(:,1)) .* (scores(:,2)>scores(:,3)).* (scores(:,2)>scores(:,4));
        elseif(c==2)
            classification = (scores(:,3)>scores(:,1)) .* (scores(:,3)>scores(:,2)).* (scores(:,3)>scores(:,4));
        else
            classification = (scores(:,4)>scores(:,1)) .* (scores(:,4)>scores(:,2)).* (scores(:,4)>scores(:,3));
        end
        
        precision(c,1)=compute_precision(classification,lab);
        accuracy(c,1)=compute_accuracy(classification,lab);
        sensitivity(c,1)=compute_sensitivity(classification,lab);
        specificity(c,1)=compute_specificity(classification,lab);
        DSC(c,1)=compute_DSC(classification,lab);
        
        %if not, assign dummy value of -1
    else
        AUC(c,1)=-1;
        precision(c,1)=-1;
        accuracy(c,1)=-1;
        sensitivity(c,1)=-1;
        specificity(c,1)=-1;
        DSC(c,1)=-1;
    end
end

return
end

function s = compute_sensitivity(classification,labels)

TP = sum(classification.*labels);
FN = sum(labels.*(1-classification));

s = TP/(TP+FN);

return
end

function s = compute_specificity(classification,labels)

TN = sum((1-classification).*(1-labels));
FP = sum(classification.*(1-labels));

s = TN/(TN+FP);

return
end

function s = compute_accuracy(classification,labels)

TN = sum((1-classification).*(1-labels));
FP = sum(classification.*(1-labels));

s = TN/(TN+FP);

return
end

function s = compute_precision(classification,labels)

TP = sum(classification.*labels);
FP = sum(classification.*(1-labels));

s = TP/(TP+FP);

return
end

function DSC = compute_DSC(classification,lab)

DSC=2*sum(classification.*lab)/(sum(classification)+sum(lab));

return
end
