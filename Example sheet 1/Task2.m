%--------------------------
%     Example Sheet 1
%         Task 2
%           a.)
%--------------------------

%--------------------------
%     Initialization
%--------------------------
clc, clear all
nIterations=5;

nFrequency=0;
nFrequencyMax=9;
probOfFixedPopulationsAtA1=zeros(nFrequencyMax,2);
probOfFixedPopulationsAtA2=zeros(nFrequencyMax,2);
frequencies=0.1:0.1:0.9;
numberOfGenerationsItTakesToFixateA1=zeros(nFrequencyMax,2);
numberOfGenerationsItTakesToFixateA2=zeros(nFrequencyMax,2);
numberOfTimesA1FixatesPopulation=zeros(nFrequencyMax,2);
numberOfTimesA2FixatesPopulation=zeros(nFrequencyMax,2);
%--------------------------
%       Simulation
%--------------------------
tic
for frequency=0.1:0.1:0.9
  nFrequency=nFrequency+1;
  
  for population=1:1
    
    if population==1
      N=100;
    elseif population==2
      N=5000;
    end
    
    
    for iteration=1:nIterations
      
      individual=zeros(N/2,2);
      individual2=zeros(N/2,2);
      
      %--------Generation of initial population--------
      
      randomNumber=rand(N/2,1);
      randomNumber2=rand(N/2,1);
      for nIndividual=1:N/2
        
        if frequency>randomNumber(nIndividual,1) % Gene 1
          individual(nIndividual,1)=1;
        elseif frequency<randomNumber(nIndividual,1)
          individual(nIndividual,1)=2;
        else
          disp('Error genotype generator allele 1')
        end
        
        if frequency>randomNumber2(nIndividual,1) %gene 2
          individual(nIndividual,2)=1;
        elseif frequency<randomNumber2(nIndividual,1)
          individual(nIndividual,2)=2;
        else
          disp('Error genotype generator allele 2')
        end
        
      end
      %----------------------------------
      %             Breeding
      %----------------------------------
      fixatedPopulation=0;
      numberOfGenerations=0;
      
      
      
      while fixatedPopulation==0
        numberOfGenerations=numberOfGenerations+1;
        
        %Selection of parents and breeding
        r=randi([1 N/2],N/2,2); % random parents
        r2=randi([1 2],N/2,2); % random genes from parents
        for i= 1:N/2
          individual2(i,:)=[individual(r(i,1),r2(i,1)) individual(r(i,2),r2(i,2))];
        end
        
        individual=individual2;
        
        if sum(all(individual<1.5))==2  % all genes are A_1
          fixatedPopulation=1;
          probOfFixedPopulationsAtA1(nFrequency,population)=probOfFixedPopulationsAtA1(nFrequency,population)+1/nIterations;
          
          
          numberOfGenerationsItTakesToFixateA1(nFrequency,population)=...
            numberOfGenerationsItTakesToFixateA1(nFrequency,population)+numberOfGenerations;
          
          numberOfTimesA1FixatesPopulation(nFrequency,population)=...
            numberOfTimesA1FixatesPopulation(nFrequency,population)+1;
          
          
          
        elseif sum(all(individual>1.5))==2   % all genes are A_2
          fixatedPopulation=1;
          
          probOfFixedPopulationsAtA2(nFrequency,population)=probOfFixedPopulationsAtA2(nFrequency,population)+1/nIterations;
          
          
          numberOfGenerationsItTakesToFixateA2(nFrequency,population)=...
            numberOfGenerationsItTakesToFixateA2(nFrequency,population)+numberOfGenerations;
          
          numberOfTimesA2FixatesPopulation(nFrequency,population)= numberOfTimesA2FixatesPopulation(nFrequency,population)+1;
        end
        
        
        
      end
      
    end
    
  end
  
end
time=toc;
timeInHours=time/nIterations*1000/60/60;
%%
%------------------------------
%            PLOT
%------------------------------

hold on
plot(frequencies,probOfFixedPopulationsAtA1(:,1),'ko--')
plot(frequencies,probOfFixedPopulationsAtA1(:,2),'ro--')
plot(frequencies,frequencies,'b')
legend('A_1,N=100','A_1,N=1000','Corresponding theory')
xlabel('Frequency')
ylabel('Probability of A_1 fixates the population')
title('Probability of fixating A1 with respect to initial frequency of A1')


%%

%--------------------------
%     Example Sheet 1
%         Task 2
%           b.)
%--------------------------
for i=1:nFrequencyMax
avgTimeToFixationOfA1(i)=...
numberOfGenerationsItTakesToFixateA1(i,1)/numberOfTimesA1FixatesPopulation(i,1);% only check for pop=100
avgTimeToAnyFixation(i)=...
  (numberOfGenerationsItTakesToFixateA1(i,1)+numberOfGenerationsItTakesToFixateA2(i,1))...
  /(numberOfTimesA1FixatesPopulation(i,1)+numberOfTimesA2FixatesPopulation(i,1));
end
%% Plot 2b
N=100;
for i=1:size(frequencies,2)
tA1(i)=-2*N*((1-frequencies(i))/frequencies(i)*log(1-frequencies(i)));
end
for i=1:size(frequencies,2)
tAny(i)=-2*N*(frequencies(i)*log(frequencies(i))+(1-frequencies(i))*log(1-frequencies(i)));
end
figure(1)
hold on
plot(frequencies,avgTimeToFixationOfA1,'ro--')
plot(frequencies,tA1,'g' )
xlabel('Frequency')
ylabel('Avrage number of generations')
legend('Simulation','Theory')
title('Avrage time for A1 to fixate population')
hold off
figure(2)
hold on
plot(frequencies,avgTimeToAnyFixation,'ro--')
plot(frequencies,tAny,'g' )
xlabel('Frequency')
ylabel('Avrage number of generations')
legend('Simulation','Theory')
title('Avrage time for any gene to fixate population')








