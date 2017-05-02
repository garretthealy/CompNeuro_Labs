%% PSET 5 
% 
% Name: Garrett Healy 
% 
% This problem set is on Principal Component Analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% In data, columns represent joint motions and rows represent 4ms time bins
% from 337 trials

% All the trials should be combined into the 22 motions
clear all 
close all 


load('Kinematics.mat')

trials = Kinematics.Trials;

% put all the data into one 22 column matrix 

data = trials{1};
for i = 2:length(trials)
    data = [data; trials{i}];
end

% getting the pca data from the trials 

s = size(data);

[coeff, score, eigs] = pca(data);
neig = flipud(sort(eigs));
sm = sum(neig);
for j = 1:length(neig)
    neig(j) = (neig(j))/sm;
end

scoeff = size(coeff);
ncoeff = zeros(scoeff(1),scoeff(2));
sscore = size(score);
nscore = zeros(sscore(1),sscore(2));

% sorting the data by eigenvalue 

for i = 1:length(neig)
    for j = 1:length(eigs)
        if neig(i)*sm == eigs(j)
            nscore(:,i) = score(:,j);
            ncoeff(:,i) = coeff(:,j);
        end
    end
end

c = 0;
p = 0;
while p<0.9 && c<22
        c = c+1;
        p = p + neig(c);
end

disp('The number of principal components it takes to explain 90% of the variance is ' + string(c) + '.')

% The first PC expalins ~50% of the variance, the second ~22%, and the
% third ~7%. 

mu = mean(score);

bounds = [min(min(score)) max(max(score))];
PCplot(coeff, mu,1,[-160 15], bounds);
PCplot(coeff, mu,2,[-160 15], bounds);
%PCplot(coeff, mu, n, [-160 15], [min(score(:, 1)) max(score(:, 1))]);

% The first principal component represents the monkey grasping at an object
% with its palm up, whereas the second represents the monkey with its paw
% down. 
%%
x = trials{1};
avg = mean(x);
reconwrflex1 = score(1:length(x),1)*coeff(:,1)'+avg;
reconwrflex2 = score(1:length(x),1:2)*coeff(:,1:2)' + avg;
reconwrflex3 = score(1:length(x),1:3)*coeff(:,1:3)' +avg;

t = 1:length(x);

figure 
hold on
plot(t,x(:,3));
plot(t,reconwrflex1(:,3));
plot(t,reconwrflex2(:,3));
plot(t,reconwrflex3(:,3));
legend('Original Data','1 PC','2 PCs','3 PCs');
hold off

% By adding more PC's, we get closer to modeling the original data. Just
% using 3 PC's gets very close. 

%%

proj1 = score(1:length(x),1);
proj2 = score(1:length(x),2);
proj3 = score(1:length(x),3);

figure
hold on
subplot(1,2,1);
hold on
plot(t,proj1);
plot(t,proj2);
plot(t,proj3);
legend('1 PC Projection','2 PC Projection','3 PC Projection');
xlabel('Observations');ylabel('Projections');
hold off
subplot(1,2,2);
plot3(proj1,proj2,proj3);
xlabel('First PC');ylabel('Second PC');zlabel('Third PC');
hold off

proj19 = score(1:length(x),19);
proj20 = score(1:length(x),20);
proj21 = score(1:length(x),21);

figure
hold on
subplot(1,2,1);
hold on
plot(t,proj19);
plot(t,proj20);
plot(t,proj21);
legend('19th PC Projection','20th PC Projection','21st PC Projection');
xlabel('Observations');ylabel('Projections');
hold off
subplot(1,2,2);
plot3(proj19,proj20,proj21);
xlabel('19th PC');ylabel('20th PC');zlabel('21st PC');
hold off

% When looking at the first few PC's, they seem to have a well defined path
% both in the observations and in the 3-d PC space. In contrast, the last
% few PC's oscillate wildly in short observation changes, indicating that
% they posisbly include a lot of noise. In addition, the amplitude of the
% lesser PC projections is extremely low, moving between -3 and 3, whereas
% the projections that account for more variance exist between values
% lesser and greater than -50 and 50. 




















%%%%%%%%%%%%%%%%%%%%%End of Script%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%