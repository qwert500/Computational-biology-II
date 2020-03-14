%-------------------------
%      Comp. Bio II
%     Example sheet 4
%         Task 1
%-------------------------
clc, clear
%-------------------------
%     Initialization
%-------------------------
N=1000;
s=[10^(-4) 5*10^(-4) 10^(-3) 5*10^(-3) 10^(-2) 5*10^(-2) 10^(-1)];
nGenes=2;
nIterations=1000;
p0=1/(2*N);

newDiploids=zeros(2*N,1);
probabilityOfFixatingAA=zeros(size(s,2),1);
%-------------------------
%      Simulation
%-------------------------
tic
for sIndex=1:size(s,2) % for all values of s
  w11=1;
  w12=1-s(sIndex)/2;
  w22=1-s(sIndex);
  for iteration=1:nIterations
    %------------------------
    % Generating population
    %------------------------
    diploid=zeros(N,nGenes);
    diploid(1,1)=1; % One geneotype is Aa
    while (sum(sum(diploid))<2*N) && (sum(sum(diploid))>0)
      
      f11=size(find(sum(diploid,2)==2),1)*(1/size(diploid,1)); % number of AA
      f12=size(find(sum(diploid,2)==1),1)*(1/size(diploid,1)); % number of Aa
      f22=size(find(sum(diploid,2)==0),1)*(1/size(diploid,1)); % number of aa
      
      wTotal=w11*f11+w12*f12+w22*f22;
      
      proportionOfAA=w11*f11/wTotal;
      proportionOfAa=w12*f12/wTotal;
      proportionOfaa=w22*f22/wTotal;
      sumOfProportions=proportionOfAA+proportionOfAa+proportionOfaa;
      
      randomNumber=sumOfProportions*rand(2*N,1);
      %--------------------------
      %          Mating
      %--------------------------
      for i=1:size(randomNumber,1)
        if randomNumber(i)<proportionOfAA
          newDiploids(i)=1;
        elseif randomNumber(i)>proportionOfAA && (randomNumber(i)<(proportionOfAA+proportionOfAa))
          newDiploids(i)=randi([0,1],1);
        else
          newDiploids(i)=0;
        end
      end
      diploid=[newDiploids(1:N) newDiploids(N+1:2*N)];
    end
    %------------------------
    %    Collecting data
    %------------------------
    if sum(sum(diploid))==2*N
      probabilityOfFixatingAA(sIndex)=probabilityOfFixatingAA(sIndex)+1/nIterations;
    end
  end
end
toc
%% PLOTS
Pkimura=(1-exp(-2*s*N*p0))./(1-exp(-2*s*N));
Phaldane=s;

hold on
plot(log(s),probabilityOfFixatingAA,'--ob')
plot(log(s),Pkimura,'m')
plot(log(s),Phaldane,'g')
xlabel('s')
ylabel('Probability for A to fixate population')
title('Probability of fixating A compared with theory from Kimura and haldane ')
legend('Simulation','Kimura','Haldane')
%%
%-------------------------
%      Comp. Bio II
%     Example sheet 4
%         Task 1
%       2:a delen
%-------------------------
clc, clear
%-------------------------
%     Initialization
%-------------------------
Ns=[2*10^2 5*10^2 10^3 5*10^3 10^4 10^5];
s=10^(-3);
nGenes=2;
nIterations=30000;


w11=1;
w12=1-s/2;
w22=1-s;
probabilityOfFixatingAA=zeros(size(Ns,2),1);
%-------------------------
%      Simulation
%-------------------------
tic
for NIndex=1:size(Ns,2) % for all values of s
    N=Ns(NIndex);
    newDiploids=zeros(2*N,1);

  for iteration=1:nIterations
    %------------------------
    % Generating population
    %------------------------
    diploid=zeros(N,nGenes);
    diploid(1,1)=1; % One geneotype is Aa
    while (sum(sum(diploid))<2*N) && (sum(sum(diploid))>0)
      
      f11=size(find(sum(diploid,2)==2),1)*(1/size(diploid,1)); % number of AA
      f12=size(find(sum(diploid,2)==1),1)*(1/size(diploid,1)); % number of Aa
      f22=size(find(sum(diploid,2)==0),1)*(1/size(diploid,1)); % number of aa
      
      wTotal=w11*f11+w12*f12+w22*f22;
      
      proportionOfAA=w11*f11/wTotal;
      proportionOfAa=w12*f12/wTotal;
      proportionOfaa=w22*f22/wTotal;
      sumOfProportions=proportionOfAA+proportionOfAa+proportionOfaa;
      
      randomNumber=sumOfProportions*rand(2*N,1);
      %--------------------------
      %          Mating
      %--------------------------
      for i=1:size(randomNumber,1)
        if randomNumber(i)<=proportionOfAA
          newDiploids(i)=1;
        elseif randomNumber(i)>proportionOfAA && (randomNumber(i)<=(proportionOfAA+proportionOfAa))
          newDiploids(i)=randi([0,1],1);
        else
          newDiploids(i)=0;
        end
      end
      diploid=[newDiploids(1:N) newDiploids(N+1:2*N)];
    end
    %------------------------
    %    Collecting data
    %------------------------
    if sum(sum(diploid))==2*N
      probabilityOfFixatingAA(NIndex)=probabilityOfFixatingAA(NIndex)+1/nIterations;
    end
  end
end
toc/60
%% PLOTS
p0=1./(2*Ns);
Pkimura=(1-exp(-2*s*Ns.*p0))./(1-exp(-2*s*Ns));
Phaldane=ones(size(Ns,2),1)*s;
hold on
plot(log(Ns),probabilityOfFixatingAA,'--ob')
plot(log(Ns),Pkimura,'m')
plot(log(Ns),Phaldane,'g')
%plot(log(Ns),log(Pkimura'+Phaldane),'k')
xlabel('log(N)')
ylabel('log(Probability of fixating AA)')
title('Fixation probability of AA, averaged over 10000 iterations')
legend('simulation','Theory from Haldane','Theory from Haldane summed with Kimura')


