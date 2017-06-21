function [avgRI,stdRI,avgRI_random,stdRI_random] = fcm_stability(pulse,ks2try)
%FCM_STABILITY Performs FCM on area response curves of input PULSEs many
%times (default = 1000) and reports the average rand_index.
%
% [RI,stdRI,randomRI,randomRI_std] = pulses.fcm_stability(ks2try)

% Process data to cluster (throw out ones > 1 NaN and turn all remaining
% NaN into 0
fits = pulse.getFits;
X = cat(1,fits.corrected_area_norm);
X = X(  sum(isnan(X),2) < 1, : );
X( isnan(X) ) = 0;
X = bsxfun(@minus,X,mean(X));
X = bsxfun(@rdivide,X,std(X));

% FCM options
o = [2 1000 1e-5 0];

for k = ks2try
    
    Niter = 1000;
    labels_all = nan( Niter, size(X,1) );
    labels_rand = nan( Niter, size(X,1) );
    
    for i = 1:Niter
        [~,U] = fcm(X,k,o);
        [~,labels_all(i,:)] = max(U);
        labels_rand(i,:) = randi(k,size(X,1),1);
        
        if mod(i,Niter/10) == 0, display('.'); end
    end
    
    RI = zeros(Niter);
    RI_random = zeros(Niter);
    
    for i = 1:Niter
        for j = 1:Niter
            RI(i,j) = valid_RandIndex( labels_all(i,:), labels_all(j,:) );
            RI_random(i,j) = valid_RandIndex( labels_rand(i,:),labels_rand(j,:) );
        end
    end
    
    display(['Done with k = ' num2str(k) ' clusters']);
    avgRI(k-1) = mean(RI(:));
    stdRI(k-1) = std(RI(:));
    
%     [s,h] = silhouette(X,labels_all(1,:));
    
%     sil{k-1} = mean(s);
    avgRI_random(k-1) = mean(RI_random(:));
    stdRI_random(k-1) = std(RI_random(:));
    
end

errorbar(2:10,avgRI,stdRI),xlabel('# of clusters'),ylabel('Rand index')
hold on,errorbar(2:10,avgRI_random,stdRI_random,'r-')

end