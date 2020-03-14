%--------------------------
%     Example Sheet 1
%         Task 3
%           b.)
%--------------------------

%--------------------------
%     Initialization
%--------------------------
clc, clear all
nIterations=1000;
nFrequency=0;
nFrequencyMax=1;
theta=[0 0.1 0.5 1 2 5];
nPopulations=size(theta,2); %given for every theta
initialFrequency=0.5;
N=1000;
u=theta./(2*N);
nGenerations=10*N;

H2=zeros(size(theta,2),nGenerations);
numberOfAllelesPerIndividual=1;
nAlleles=2;
%--------------------------
%       Simulation
%--------------------------
tic


for population=1:nPopulations
  for iteration=1:nIterations
    
    individual=zeros(N,numberOfAllelesPerIndividual);
    individual2=zeros(N,numberOfAllelesPerIndividual);
    
    %--------Generation of initial population--------
    
    randomNumber=rand(N,numberOfAllelesPerIndividual);
    randomNumber2=rand(N,numberOfAllelesPerIndividual);
    
    for nIndividual=1:size(individual,1) %generates genes for individuals
      for gene=1:size(individual, 2)
        if initialFrequency>randomNumber(nIndividual,gene)
          individual(nIndividual,gene)=1;
        elseif initialFrequency<randomNumber(nIndividual,gene)
          individual(nIndividual,gene)=2;
        else
          disp('Error genotype generator in initial step')
        end
      end
    end
    %----------------------------------
    %             Breeding
    %----------------------------------
    for generation=1:nGenerations
      r=randi([1 N],N,1);
      r2=rand(N,1);
      individual2=individual(r(:));

      mutationIndex=find(r2<u(population));
      if population==1
      else
      for j=1:size(mutationIndex,1)
        if individual2(mutationIndex(j))==1
          individual2(mutationIndex(j))=2;
        elseif individual2(mutationIndex(j))==2
          individual2(mutationIndex(j))=1;
        else
          print('ERROR in mutation')
        end
      end
      end
      individual=individual2;
      nA1=length(find ( individual == 1 ));
      nA2=N-nA1;
      if nA1<2 || nA2<2
        H2(population,generation)=H2(population,generation)+0;
      else
      H2(population,generation)=H2(population,generation)+...
        (1-((nchoosek(nA1,nAlleles)+nchoosek(nA2,nAlleles))/nchoosek(nA1+nA2,nAlleles)))/nIterations;
      end 
    end
  end
end
time=toc;
timeInHours=time/nIterations*1000/60/60;
load gong 
sound(y,Fs)
%%
%-------------------------
%         Plots           
%-------------------------
H2theory=zeros(size(theta,2),nGenerations);
for i=1:size(theta,2)
H2theory(i,:)=theta(i)/(1+2*theta(i));
end


meanSteadyState=mean(H2,2);
meanSteadyStateTheory=(theta)./(1+2*theta);
figure(1)
for k=1:size(H2,1)
subplot(3,2,k)
hold on
plot(1:1:nGenerations,H2(k,:),'b')
plot(1:1:nGenerations,H2theory(k,:),'g')
xlabel('Generations')
ylabel('Hetrozygosity')
title(['Probability of hetrozygosity,','u=',num2str(u(k))])
legend('simulation','Theory')
end
hold off
figure(2)
hold on 
plot(theta,meanSteadyState,'ro--')
plot(theta,meanSteadyStateTheory,'g')
xlabel('\theta')
ylabel('Mean of Hetrozygosity in steady state')
title('Mean of hetrozygosity at steady state depending on the parameter \theta')
legend('Simulation', 'Theory')
%%
%--------------------
%      Task 3 c.)
%--------------------
clc
for i=2:4
  subplot(2,2,i)
histogram(H2(i,2000:nGenerations))
xlabel('Value of hetrozygosity')
ylabel('number of data points ')
title(['Histogram of hetrozygosity','\theta=',num2str(theta(i))])
end
for i=2:4
variance(i-1)=var(H2(i,2000:nGenerations)); % at SS
end
for i=2:4
Mean(i-1)=mean(H2(i,2000:nGenerations));% at SS
end
