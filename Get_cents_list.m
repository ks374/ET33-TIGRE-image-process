function [cents] = Get_cents_list(stats)
%get the weighted centroid list from the structure stats
    cents = zeros(numel(stats),3);
    for i = 1:numel(stats)
        cents(i,:) = stats(i).WeightedCentroid;
    end
end

