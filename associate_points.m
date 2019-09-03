function [locs1_, locs2_] = associate_points(locs1, locs2, thresh)
%% Associate peaks 

if ~exist('thresh', 'var') || isempty(thresh)
    thresh = 0.1;
end

% use pdist2 but break into windows
fprintf('associating peaks...');
dloc = median(diff(locs1));
locs1_ = [];
locs2_ = [];
N = min(length(locs1), 200);  % This is the window size for computational efficiency
for ii = 1:N:length(locs1)-N+1    
    inds = locs2 >= locs1(ii) - 2*dloc & locs2 <= locs1(ii + N - 1) + 2*dloc;
    if nnz(inds) == 0
        continue
    end
    locs1b = locs1(ii:ii+N-1);
    locs2b = locs2(inds);    
    d = pdist2(locs1b, locs2b);    
    dmatch = d < thresh * dloc;
    % assert(max(sum(dmatch, 1)) <= 1 && max(sum(dmatch, 2)) <= 1);    
    
    % I need an elegant way of removing multiple pixels on the same row or column
    temp = NaN(size(d));
    temp(dmatch) = d(dmatch);
    for jj = 1:size(dmatch, 1)
        if nnz(dmatch(jj, :)) > 1
            dmatch(jj, :) = false;
            [~, mi] = min(temp(jj, :));
            dmatch(jj, mi) = true;
        end
    end
    temp = NaN(size(d));
    temp(dmatch) = d(dmatch);
    for jj = 1:size(dmatch, 2)
        if nnz(dmatch(:, jj)) > 1
            dmatch(:, jj) = false;
            [~, mi] = min(temp(:, jj));
            dmatch(mi, jj) = true;
        end
    end
    
    inds1 = any(dmatch, 2);    
    inds2 = any(dmatch, 1);
    
    locs1_ = [locs1_; locs1b(inds1)];
    locs2_ = [locs2_; locs2b(inds2)];
end

[locs1_, ia] = unique(locs1_);
locs2_ = locs2_(ia);

[locs2_, ia] = unique(locs2_);
locs1_ = locs1_(ia);

fprintf('done\n');