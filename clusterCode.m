spikes = [];
for i = 1:length(cluster_class)
    if cluster_class(i) == 1
        spikes = [spikes; cluster_class(i,2)];
    end
end
