%-------------------------
%      Computational
%       Biology II
%         Task 2a.)
%-------------------------
clc, clear
%-------------------------
%     Initialization
%-------------------------
n=2;
Rs=[0 10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 1 2 3 5 10];
nIterations=1000;

nGenes=2;
numberOfLinages=n;
diploid=ones(n,2);
tau=zeros(2,nIterations,size(Rs,2));
timeIntoThePast=0;
%-------------------------
%     algorithm
%-------------------------
tic
for nR=1:size(Rs,2)
  R=Rs(nR);
for iteration=1:1:nIterations
  numberOfLinages=n;
  timeIntoThePast=0;
  diploid=ones(n,2);
while tau(1, iteration,nR)==0 || tau(2,iteration,nR)==0
  
  %-------------------------
  %   Generating Times
  % for coalescence and
  %    recombination
  %-------------------------
  
  k=numberOfLinages;
  l=Calculatingl(diploid,n);
  
  lambdaCoalescense=nchoosek(k,2);
  lambdaRecombination=l*R/2;
  
  TCoalescense=exprnd(1/lambdaCoalescense);
  TRecombination=exprnd(1/lambdaRecombination);
  %-------------------------
  %     Recombination or
  %       Coalescense
  %-------------------------
  
  
  if TRecombination<TCoalescense ...
      && (size(find(diploid==0),1)+size(find(diploid==n),1))<size(diploid,1)
    condition=0;
    while condition==0 % find diploid that are able to recombine
      recombination=randi([1 size(diploid,1)],1,1);
      
      if size(find(diploid(recombination,:)==0),2)==0 ...
          && size(find(diploid(recombination,:)==n),2)==0
        
        condition=1;
      end
    end
    
    newDiploid1=[diploid(recombination,1) 0];
    newDiploid2=[0 diploid(recombination,2)];
    
    diploid(recombination,:)=[];
    diploid=[diploid;newDiploid1;newDiploid2];
    timeIntoThePast=TRecombination+timeIntoThePast;
    
  else %THE COALESCENSE PROCESS
    coalescense=coalescenseOrphans(diploid); % random coalescense
    newDiploid=diploid(coalescense(1),:)+diploid(coalescense(2),:);
    
    if coalescense(1)<coalescense(2)% killing the correct children
      diploid(coalescense(1),:)=[];
      diploid(coalescense(2)-1,:)=[];
    else
      diploid(coalescense(1),:)=[];
      diploid(coalescense(2),:)=[];
    end
    
    diploid=[diploid; newDiploid]; % Saving Ancestor
    timeIntoThePast=TCoalescense+timeIntoThePast;
  end
  
  % for i=1:1:size(diploid,1) % sets gene to the correct number,i.e to 2 if >2
  %   for j=1:1:nGenes
  %     if diploid(i,j)>n
  %       diploid(i,j)=n;
  %     end
  %   end
  % end
  %-------------------------
  %     Time to MRCA
  %-------------------------
  
  for i=1:size(diploid,1) % Check if MRCA is found for any gene and sets the time
    for j=1:nGenes
      if diploid(i,j)==n && tau(j,iteration,nR)==0
        
        tau(j,iteration,nR)=timeIntoThePast;
      end
    end
  end
  
  numberOfLinages=size(diploid,1);
end
end
end
toc
%% Plots and histograms
tauA=zeros(nIterations,size(Rs,2));
tauB=zeros(nIterations,size(Rs,2));
for nR=1:1:size(Rs,2)
tauA(:,nR)=tau(1,:,nR);
tauB(:,nR)=tau(2,:,nR);

end
x=linspace(0,10);
% choose nR for every plot 
subplot(1,2,1)
hold on
histogram(tauA(:,6),50)
plot(x,exppdf(x,1)*nIterations/5,'g','Linewidth',2)
ylabel('Frequency')
title('Histogram for \tau_a plotted with theoretical distributaion, R=0.1')
subplot(1,2,2)
hold on
histogram(tauB(:,6),50)
plot(x,exppdf(x,1)*nIterations/5,'g','Linewidth',2)
ylabel('Frequency')
title('Histogram for \tau_b plotted with theoretical distributaion, R=0.1')
%-----------------------
% Means from Simulation
%-----------------------
for nR=1:size(Rs,2)
meanOfTauA(nR)=mean(tauA(:,nR));
meanOfTauB(nR)=mean(tauB(:,nR));
meanOfTauAsquare(nR)=mean(tauA(:,nR).^2);
meanOfTauBsquare(nR)=mean(tauB(:,nR).^2);
end
%% 
%-------------------------
%      Computational
%       Biology II
%         Task 2b.)
%-------------------------      
% Use values from task a!
%-------------------------
%    Initialization
%-------------------------
tauA=zeros(nIterations,size(Rs,2));
tauB=zeros(nIterations,size(Rs,2));

for nR=1:1:size(Rs,2)
tauA(:,nR)=tau(1,:,nR);
tauB(:,nR)=tau(2,:,nR);
TAU(nR)=mean(tauA(:,nR))./2+mean(tauB(:,nR))./2;
TAUsquare(nR)=mean(tauA(:,nR).^2)./2+mean(tauB(:,nR).^2)./2;
end
for nR=1:1:size(Rs,2)
  correlationFunction(nR)=(mean(tauA(:,nR).*tauB(:,nR))-mean(tauA(:,nR))*mean(tauB(:,nR)))/...
    (TAUsquare(nR)-TAU(nR)^2);
  theoreticalCorrelationFunction(nR)=(Rs(nR)+18)/(Rs(nR)^2+13*Rs(nR)+18);
end

hold on

plot(log(Rs),correlationFunction,'bo--')
plot(log(Rs),theoreticalCorrelationFunction,'g')
xlabel('log(R)')
ylabel('Correlation function')
title('Correlation function as a function of R')
legend('Simulation','Theory')












