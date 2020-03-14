%--------------------------
%       Computational
%        Biology II
%     Examples sheet 2
%          Task 1
%   RoR INTE DENNA MER NU!
%--------------------------
clc, clear
%--------------------------
%      Initialization
%--------------------------
%n=2; % sample size task 2
%theta=5; % sample size task 2
n=100; % sample size task 1
theta=[1 10 100]; % sample task 1
N=1000;
nIterations=1000;

numberOfDiffrentAlleles=zeros(size(theta,2),nIterations);
u=theta./(2*N);
lambda=zeros(n-1,1);
T=zeros(n-1,1);
genealogy=zeros(n+n-1,7);
numberOfExistingNodes=n;
listOfOrphans=(1:1:n)';
timeIntoThePast=0;
totalNumberOfMutations=zeros(nIterations,1);

for iteration=1:nIterations
  for nTheta=1:size(theta,2)

    
    for i=2:n
      lambda(i)=nchoosek(i,2)/N;
    end
    T=exprnd((lambda(2:end)).^(-1));
    T=sort(T);% smallest value first ...
    %--------------------------------
    %     Generating MRCA and
    %          Dendrogram
    %--------------------------------
    tic
    for i=1:1:n-1 % generating MRCA from n samples from population
      if size(listOfOrphans,1)>1
        coalescenseOrphan=coalescenseOrphans(listOfOrphans);
        
        genealogy(listOfOrphans(coalescenseOrphan(1)),3)=numberOfExistingNodes+1; % parents
        genealogy(listOfOrphans(coalescenseOrphan(2)),3)=numberOfExistingNodes+1;
        
        genealogy(numberOfExistingNodes+1,1)=listOfOrphans(coalescenseOrphan(1));% Descendants
        genealogy(numberOfExistingNodes+1,2)=listOfOrphans(coalescenseOrphan(2));
        
        timeIntoThePast=T(i)+timeIntoThePast;
        genealogy(numberOfExistingNodes+1,4)=timeIntoThePast; % Time for coalescent event
        
        numberOfExistingNodes=numberOfExistingNodes+1;
        
        listOfOrphans(coalescenseOrphan(1))=[];
        if coalescenseOrphan(1)<coalescenseOrphan(2)
          listOfOrphans(coalescenseOrphan(2)-1)=[];
        else
          listOfOrphans(coalescenseOrphan(2))=[];
        end
        listOfOrphans=[listOfOrphans;numberOfExistingNodes];
      else
      end
    end
    
    %----------------------
    %      Mutations
    %----------------------
    T=sort(T,'descend');
    
    for i=n+1:n+n-1
      for j=1:2
        timeToDiscard=genealogy(genealogy(i,j),4);
        genealogy(i,4+j)=poissrnd(u(nTheta)*(genealogy(i,4)-timeToDiscard)); % time on each branch
      end
    end
    
    %----------------------
    %   Counting diffrent
    %       alleles
    %----------------------
    nNewGenes=0;
    for i=0:n+n-2
      for j=1:2
        if genealogy(n+n-1-i,4+j)>0
          nNewGenes=nNewGenes+1;
          genealogy(genealogy(n+n-1-i,j),7)=nNewGenes;
        elseif genealogy(n+n-1-i,j)==0
          % No more Descendants
        else
          genealogy(genealogy(n+n-1-i,j),7)=genealogy(n+n-1-i,7);
        end
      end
    end
    numberOfDiffrentAlleles(nTheta,iteration)=size(unique(genealogy(:,7)),1);
  end
totalNumberOfMutations(iteration)=sum(genealogy(:,5))+sum(genealogy(:,6));
end

%%
%-------- Theory & Plots Task 1 --------%
load stirlingsNumber.mat %Taken from Mathematica
stirlingNumber=stirlingsNumber;
m=1:n; %number of possible allelic types

%Ewens sampling formula
nominator=zeros(length(theta),n);
denominator=zeros(length(theta),n-1);
P_n=zeros(size(theta,2),n);
K=zeros(length(theta),1);

for p=1:size(theta,2) %for all different theta
    for k=2:n
        denominator(p,k-1)=(k-1)+theta(p);
    end

    nominator(p,:)=abs(stirlingsNumber).*theta(p).^(m-1);    
    P_n(p,:)=nominator(p,:)./prod(denominator(p,:),2);
end

for i=1:size(theta,2)
    hold on
  subplot(1,3,i)
      hold on
    histogram(numberOfDiffrentAlleles(i,:),'Facecolor','yellow')
  plot(1:1:n,P_n(i,:)*nIterations,'g','Linewidth',2)
  xlabel('Number of allelic types')
  ylabel('Frequency of observing m allelic types')
  title('Histogram  with theory, avraged over 1000 iterations.')
  legend('Simulation',['Theory, \theta=',num2str(theta(i))])
end



