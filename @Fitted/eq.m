function equality = eq(fit1,fit2)
% Equality comparator for FITTED
% right now slow, expecting array in first argument and single object
% in the second (fit2). Will return equal if the width_frame of two
% fits have overlap > 3.
if numel(fit1) > 1 && numel(fit2) > 1
    error('Cannot handle TWO array inputs.');
end
%             names = setdiff(fieldnames(fit2),{'fitID','category'});
equality = false(1,numel(fit1));
for j = 1:numel(fit1)
    if fit1(j).stackID == fit2.stackID
        % can't use bsxfun because of un-uniform output
        if numel(fit1(j).width_frames( ...
                ismember(fit1(j).width_frames, fit2.width_frames))) > 5
            equality(j) = 1;
        end
    end
end

end %eq