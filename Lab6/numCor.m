function bool = isCor(i,j,k)


out = zeros(4,1); %4 bandwiths
for g = 1:4 
    neur2 = spikes.spikes{g}; %neuron to test against
    for h = 1:4
        band2 = neur2{h}; %""
        for f = 1:5
            if g~=i||h~=j||f~=k
                amp2 = band2{f}; %""
                out(h) = out(h) + spkd(amp, amp2,q);
            end
        end
        if h == j && g == i
            out(h) = out(h)/4;
        else 
            out(h)= out(h)/5;
        end
    end
end
out = out./(16);
[m,loc] = min(out);
if loc == j
    pcorrect = pcorrect+1;
end

