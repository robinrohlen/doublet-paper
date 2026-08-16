function locs_new=interpolate_spike_train(locs,source)

y=interp1(linspace(0,size(source,2)/2e3,size(source,2)),source,linspace(0,size(source,2)/2e3,20*size(source,2)),'spline');

locs_new=zeros(size(locs));

for ind=1:length(locs_new)
    [~,maxVal]=max(y(20*locs(ind)+(-20:20)));
    locs_new(ind)=20*locs(ind)+maxVal;
end