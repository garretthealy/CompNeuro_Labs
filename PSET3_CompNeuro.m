%% Problem Set # 3 
% 
% Name: Garrett Healy 
% 
% This week, we will be looking at retinal data from Michael Berry’s lab at 
% Princeton University. In this experiment, a salamander retina was isolated 
% and placed on an array of electrodes, from which the responses of 115 
% retinal ganglion cells were simultaneously recorded. White noise was 
% delivered to the retina for about an hour. We have provided you with the 
% responses from 10 neurons. You will use this white noise stimulus and the 
% evoked responses to compute and model receptive fields for these ten neurons.
% 
%% 1
% Plot a raster of the first minute of responses from each neuron in a 
% single plot (each row corresponding to a different neuron).

% 1, 4, 11, 15, 26, 51, 80, 84, 96, 105

load('retinaData.mat')

% First step is to get data from only those 10 neurons 
% Data * 100 = microseconds 

data = zeros(10,82013);
lengths = zeros(10,1);
c=1;

for i = 1:length(retinaData.spikes)
    if i == 1 || i == 4 || i == 11 || i == 15 || i == 26 || i == 51 || i == 80 ...
            || i == 84 || i == 96 || i == 105
        d = retinaData.spikes{i};
        for j = 1:length(d)
            data(c,j) = d(j);
            lengths(c) = length(d); %Get the lengths of the neurons data
        end
        c=c+1;
    end
end 

onemin = zeros(10,3300);
s = size(data);
oneminlengths = zeros(10,1);

for i = 1:s(1)
    c=1;
    while data(i,c)< 60e4 && c<3300
        onemin(i,c) = data(i,c) / 10000;
        c = c+1;
    end
    oneminlengths(i) = c;
end 

figure 
xlabel('Time (s)');ylabel('Neuron');
title('Spiking of First Minute Responses of Retinal Ganglion Cells');
hold on

for i = 1:10
    trial = onemin(i,:);
    for q = 1:oneminlengths(i)-1
        spkx = [trial(q), trial(q)];
        spky = [i - 0.4,i + 0.4];
        line(spkx,spky,'color','k','LineWidth',1);
    end
end
hold off

%% 2A
% 
% Use the white noise stimulus and the spikes from neuron #1 to write a 
% function that computes a spatiotemporal receptive field. The stimulus is 
% a 40 by 40 grid in which each checkerboard spot is either on or off. The 
% checkerboard pattern changes at a rate of 30 Hz. Your STRF should range 
% cover lags between 300 ms and -10 ms, with a resolution of 10 ms. Plot 
% your STRF on a 6 by 6 grid of subplots.
% 
% Only using neuron #1, deltastimulus = 3.34 ms, starts at 51243,
% ends at 40292744
%

close all; 

ft = retinaData.stimulusFrameTimes; % time 
sf = retinaData.stimulusFrames; % spike or no spike
one = data(1,:); % first neuron only 
s = size(sf); 
l3 = round(length(ft)/3); % 10ms resolution
counts = zeros(l3,1); %records counts for all time with 10ms res

for i = 2:length(ft)/3
    count = 0;
    high = ft(i*3);
    low = ft((i-1)*3);
    if low ~=0
        for j = 1:lengths(1) %if there is a spike in the 10ms time window, this counts it
            if one(j)>low && one(j)<=high
                count = count+1;
            else
                counts(i-1) = count;
            end
        end 
    end
end 

heatvals = zeros(1600,length(counts));
m = max(counts);
for i = 1:length(counts)
    c = 1;
    for j = 1:s(1)
        for k = 1:s(2)
            if sf(j,k,i) == 1 
                heatvals(c,i) = counts(i)/(m);
            elseif sf(j,k,i) == 0  
                heatvals(c,i) = - counts(i)/(m);
            end       
            c = c+1;
        end
    end 
end

figure 
d = zeros(40,40);
for j = 1:32
    for i = 1:1600
        xval = ceil(i/40);
        yval = i - (40*(xval-1));
        d(xval,yval) = heatvals(i,j);
    end
    hold on
    subplot(6,6,j)
    hold on
    for i = 1:40
        for k = 1:40
            point = d(i,k);
            if point>0 
                ppoint = point;
                npoint = 0;
            elseif point<0
                ppoint = 0;
                npoint = -point;
            else 
                ppoint = 0;
                npoint = 0;
            end
            x = [i,i];
            y = [k-0.5,k+0.5];
            line(x,y,'color',[(ppoint) (1-abs(point)) (npoint)],'LineWidth',7.5);
        end
    end
    hold off
    title(['STRF for Neuron 1 at '  num2str(j*10 - 20) ...
         ' ms.']);
    hold off
end
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,'Figure1.jpg')
close all;

%% 2B
% Once you have successfully written a function that can compute the STRF 
% for neuron #1, compute and save the STRF for all 10 neurons.
% 

startt = datestr(now);

s = size(sf);
l3 = round(length(ft)/3);

for p = 1:10
    trial = data(p,:);
    counts = zeros(l3,1);
    for i = 2:length(ft)/3
        count = 0;
        high = ft(i*3);
        low = ft((i-1)*3);
        if low ~=0
            for j = 1:lengths(p)
                if trial(j)>low && trial(j)<=high
                    count = count+1;
                else
                    counts(i-1) = count;
                end
            end 
        end
        if i == round(length(ft)/6)
            disp(['Halfway Point for Neuron #'  num2str(p)])
            disp(' ')
        end
    end 

    heatvals = zeros(1600,length(counts));
    m = max(counts);
    for i = 1:length(counts)
        c = 1;
        for j = 1:s(1)
            for k = 1:s(2)
                if sf(j,k,i) == 1 
                    heatvals(c,i) = counts(i)/(m);
                elseif sf(j,k,i) == 0  
                    heatvals(c,i) = - counts(i)/(m);
                end       
                c = c+1;
            end
        end 
    end

    figure 
    d = zeros(40,40);
    for j = 1:32
        for i = 1:1600
            xval = ceil(i/40);
            yval = i - (40*(xval-1));
            d(xval,yval) = heatvals(i,j);
        end
        hold on
        subplot(6,6,j)
        hold on
        for i = 1:40
            for k = 1:40
                point = d(i,k);
                if point>0 
                    ppoint = point;
                    npoint = 0;
                elseif point<0
                    ppoint = 0;
                    npoint = -point;
                else 
                    ppoint = 0;
                    npoint = 0;
                end
                x = [i,i];
                y = [k-0.5,k+0.5];
                line(x,y,'color',[(ppoint) (1-abs(point)) (npoint)],'LineWidth',7.5);
            end
        end
        hold off
        title(['STRF for Neuron ' num2str(p) ' at '  num2str(j*10 - 20) ...
                 ' ms.']);
        hold off
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,['Figure' num2str(p) '.jpg']);
    close all;
end

endt = datestr(now);
tchange = endt - startt;
minutes = tchange(16)*10 + tchange(17);
hours = tchange(13)*10 + tchange(14);
disp([num2str(hours)  ' Hours and '  num2str(minutes) ' Minutes' ... 
   ' have elapsed.'])

%% 3
% Using lsqcurvefit, fit each of your neuron STRFs using a time-varying 2D 
% Gaussian model. You should fit the Gaussian center, the Gaussian width and 
% 32 different amplitudes (one for each time point).
% 
% For each neuron, use a series of subplots on a single figure to plot: 1) 
% the empirically determined STRF at time lag 80 ms, 2) the receptive field
% model at time lag 80 ms, and 3) the amplitude time course from 300 ms to 
% -10 ms. Did the Gaussian successfully fit every neuron?
% 
% What is the average receptive field width? What percentage of the neurons 
% in this sample show OFF response properties? What is the average time to
% peak amplitude for OFF cells? What is the average time to peak amplitude 
% for ON cells?
% 

I=imread('Figure1.jpg');%assume gray scale, not RGB
[n,m]=size(I);%assumes that I is a nxm matrix
[X,Y]=meshgrid(1:n,1:m);%your x-y coordinates
x(:,1)=X(:); % x= first column
x(:,2)=Y(:); % y= second column
f=I(:); % your data f(x,y) (in column vector)
%--- now define the function in terms of x
%--- where you use x(:,1)=X and x(:,2)=Y
fun = @(c,x) c(1)+c(2)*exp(-((x(:,1)-c(3))/c(4)).^2-((x(:,2)-c(5))/c(6)).^2)
%--- now solve with lsqcurvefit
options=optimset('TolX',1e-6);
c0=[1 1 1 1 1 1]%start-guess here
cc=lsqcurvefit(fun,c0,x,f,[],[],options)
Ifit=fun(cc,x); %your fitted gaussian in vector
Ifit=reshape(Ifit,[n m]);%gaussian reshaped as matrix
% plot with surf(X,Y,Ifit) ...








