function [bool,loc,diff]= isCor(i,j,k,q)

spikes = load('spikes.mat');

bool = 0;
diff = zeros(4,4,5); %4 bandwiths
amp = spikes.spikes{i}{j}{k};
for g = 1:4 
    neur2 = spikes.spikes{g}; %neuron to test against
    for h = 1:4
        band2 = neur2{h}; %""
        for f = 1:5
            if g~=i||h~=j||f~=k
                amp2 = band2{f}; %makes sure it isnt the same spike train
                diff(g,h,f) = spkd(amp, amp2,q);
            end
        end
    end
end
diff = diff./(79);
diff = mean(diff,3);
diff = mean(diff,1);
[~,loc] = min(diff);
if loc == j
    bool = 1;
end
end

