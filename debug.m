for i=1:length(ReverseCurrent_slct)
    for j =1:length(ReverseCurrent_slct)
        if length(ReverseCurrent_slct{i}) == length(ReverseCurrent_slct{j})
%             fprintf('%d,%d\n',i,j)
            
            if sum(ReverseCurrent_slct{i} - ReverseCurrent_slct{j}) == 0
                fprintf('found:%d,%d\n',i,j)
            end
        end
    end
end

for i=1:length(ForwardCurrent_slct)
    for j =1:length(ForwardCurrent_slct)
        if length(ForwardCurrent_slct{i}) == length(ForwardCurrent_slct{j})
%             fprintf('%d,%d\n',i,j)
            
            if sum(ForwardCurrent_slct{i} - ForwardCurrent_slct{j}) == 0
                fprintf('found:%d,%d\n',i,j)
            end
        end
    end
end