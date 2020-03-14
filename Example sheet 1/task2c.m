%--------------------------
%     Example Sheet 1
%         Task 2
%           a.)
%--------------------------

%--------------------------
%     Initialization
%--------------------------
clc, clear all
nIterations=1000;
nPopulations=6;
nFrequency=0;
nFrequencyMax=1;

tA1=zeros(nPopulations,1);
tAny=zeros(nPopulations,1);
probOfFixedPopulationsAtA1=zeros(nFrequencyMax,nPopulations);
probOfFixedPopulationsAtA2=zeros(nFrequencyMax,nPopulations);
frequencies=0.1:0.1:0.9;
numberOfGenerationsItTakesToFixateA1=zeros(nFrequencyMax,nPopulations);
numberOfGenerationsItTakesToFixateA2=zeros(nFrequencyMax,nPopulations);
numberOfTimesA1FixatesPopulation=zeros(nFrequencyMax,nPopulations);
numberOfTimesA2FixatesPopulation=zeros(nFrequencyMax,nPopulations);
numberOfAllelesPerIndividual=1;
%--------------------------
%       Simulation
%--------------------------
tic
for frequency=0.5
  nFrequency=nFrequency+1;
  
  for population=1:nPopulations
    
    if population==1
      N=10;
    elseif population==2
      N=50;
    elseif population==3
      N=100;
    elseif population==4
      N=500;
    elseif population==5
      N=1000;
    elseif population==6
      N=5000;
    end
    
    
    for iteration=1:nIterations
      
      individual=zeros(N,numberOfAllelesPerIndividual);
      individual2=zeros(N,numberOfAllelesPerIndividual);
      
      %--------Generation of initial population--------
      
      randomNumber=rand(N,numberOfAllelesPerIndividual);
      randomNumber2=rand(N,numberOfAllelesPerIndividual);
      
      for nIndividual=1:size(individual,1) %generates genes for individuals
        for gene=1:size(individual, 2)
          if frequency>randomNumber(nIndividual,gene)
            individual(nIndividual,gene)=1;
          elseif frequency<randomNumber(nIndividual,gene)
            individual(nIndividual,gene)=2;
          else
            disp('Error genotype generator in initial step')
          end
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
        r=randi([1 N],N,numberOfAllelesPerIndividual); % random parents
        
        if numberOfAllelesPerIndividual==1
          individual2(:)=individual(r(:));
          
          %         elseif numberOfAllelesPerIndividual==2
          %           individual2(:)=individual(r(:));
          %
          %           r2=randi([1 numberOfAllelesPerIndividual],N,numberOfAllelesPerIndividual); % random genes from parents
          %           for i= 1:N
          %             individual2(i,:)=[individual(r(i,1),r2(i,1)) individual(r(i,2),r2(i,2))];
          %           end
        else
          print('Error in breeding')
        end
        
        individual=individual2;
        
        if sum(all(individual<1.5))==numberOfAllelesPerIndividual  % all genes are A_1
          fixatedPopulation=1;
          probOfFixedPopulationsAtA1(nFrequency,population)=probOfFixedPopulationsAtA1(nFrequency,population)+1/nIterations;
          
          
          numberOfGenerationsItTakesToFixateA1(nFrequency,population)=...
            numberOfGenerationsItTakesToFixateA1(nFrequency,population)+numberOfGenerations;
          
          numberOfTimesA1FixatesPopulation(nFrequency,population)= numberOfTimesA1FixatesPopulation(nFrequency,population)+1;
          
          
          
        elseif sum(all(individual>1.5))==numberOfAllelesPerIndividual   % all genes are A_2
          fixatedPopulation=1;
          
          probOfFixedPopulationsAtA2(nFrequency,population)=probOfFixedPopulationsAtA2(nFrequency,population)+1/nIterations;
          
          
          numberOfGenerationsItTakesToFixateA2(nFrequency,population)=...
            numberOfGenerationsItTakesToFixateA2(nFrequency,population)+numberOfGenerations;
          
          numberOfTimesA2FixatesPopulation(nFrequency,population)= numberOfTimesA2FixatesPopulation(nFrequency,population)+1;
        end
        
      end
    end
    %Theoretical expactation
    tA1(population)=-2*N*((1-frequency)/frequency*log(1-frequency));
    tAny(population)=-2*N*(frequency*log(frequency)+(1-frequency)*log(1-frequency));
    
    
  end
end
time=toc;
timeInHours=time/nIterations*1000/60/60;
%%

%--------------------------
%     Example Sheet 1
%         Task 2
%           b.)
%--------------------------
avgTimeToFixationOfA1=zeros(nFrequencyMax,nPopulations);
avgTimeToAnyFixation=zeros(nFrequencyMax,nPopulations);

for i=1:nPopulations
  avgTimeToFixationOfA1(i)=...
    numberOfGenerationsItTakesToFixateA1(1,i)/numberOfTimesA1FixatesPopulation(1,i);% only check for pop=100
  avgTimeToAnyFixation(i)=...
    (numberOfGenerationsItTakesToFixateA1(1,i)+numberOfGenerationsItTakesToFixateA2(1,i))...
    /(numberOfTimesA1FixatesPopulation(1,i)+numberOfTimesA2FixatesPopulation(1,i));
end
figure(1)
hold on
plot([10 50 100 500 1000 5000],avgTimeToFixationOfA1,'ro--')
plot([10 50 100 500 1000 5000],tA1,'g')
xlabel('Population size')
ylabel('Avrage number of generations')
title('Avrage time for A1 to fixate population')
legend('simulation','Theory')
hold off
figure(2)
hold on
plot([10 50 100 500 1000 5000],avgTimeToAnyFixation,'ro--')
plot([10 50 100 500 1000 5000],tAny,'g')
xlabel('Population size')
ylabel('Avrage number of generations')
title('Avrage time for any gene to fixate population')
legend('simulation','Theory')
hold off


