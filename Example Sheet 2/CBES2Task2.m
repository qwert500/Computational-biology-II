%--------------------------
%       Computational
%        Biology II
%     Examples sheet 2
%          Task 2
%--------------------------
clc, clear
%--------------------------
%      Initialization
%--------------------------
n=[2 100]; % sample size
theta=5;
N=1000;
nIterations=1000;

totalNumberOfMutations=zeros(size(n,2),nIterations);

for sample=1:size(n,2)

numberOfDiffrentAlleles=zeros(size(theta,2),nIterations);
u=theta./(2*N);
lambda=zeros(n(sample)-1,1);
T=zeros(n(sample)-1,1);
genealogy=zeros(n(sample)+n(sample)-1,7);
numberOfExistingNodes=n(sample);
listOfOrphans=(1:1:n(sample))';
timeIntoThePast=0; 
  
for iteration=1:nIterations
  for nTheta=1:size(theta,2)
    for i=2:n(sample)
      lambda(i)=nchoosek(i,2)/N;
    end
    T=exprnd(1./(lambda(2:end)));
    T=sort(T);% smallest value first ...

    %--------------------------------
    %     Generating MRCA and
    %          Dendrogram
    %--------------------------------
    tic
    for i=1:1:n(sample)-1 % generating MRCA from n samples from population
      
      if size(listOfOrphans,1)>1
        
        coalescenseOrphan=coalescenseOrphans(listOfOrphans);
        
        genealogy(listOfOrphans(coalescenseOrphan(1)),3)=...
          numberOfExistingNodes+1; % parents
        genealogy(listOfOrphans(coalescenseOrphan(2)),3)=...
          numberOfExistingNodes+1;
        
        genealogy(numberOfExistingNodes+1,1)=...
          listOfOrphans(coalescenseOrphan(1));% Descendants
        genealogy(numberOfExistingNodes+1,2)=...
          listOfOrphans(coalescenseOrphan(2));
        
        timeIntoThePast=T(i)+timeIntoThePast;
        genealogy(numberOfExistingNodes+1,4)=timeIntoThePast; 
        
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
    for i=n(sample)+1:n(sample)+n(sample)-1
      for j=1:2
        timeToDiscard=genealogy(genealogy(i,j),4);
        beta=u(nTheta)*(genealogy(i,4)-timeToDiscard);
        genealogy(i,4+j)=poissrnd(beta); % time on each branch
      end
    end
    
    %----------------------
    %   Counting diffrent
    %       alleles
    %----------------------
    nNewGenes=0;
    for i=0:n(sample)+n(sample)-2
      for j=1:2
        nNewGenes=nNewGenes+1;
        if genealogy(n(sample)+n(sample)-1-i,4+j)>0
          %nNewGenes=nNewGenes+1;
          genealogy(genealogy(n(sample)+n(sample)-1-i,j),7)=nNewGenes;
        elseif genealogy(n(sample)+n(sample)-1-i,j)==0
          % No more Descendants
        else
          genealogy(genealogy(n(sample)+n(sample)-1-i,j),7)=...
            genealogy(n(sample)+n(sample)-1-i,7);
        end
      end
    end
    
    numberOfDiffrentAlleles(nTheta,iteration)=...
      size(unique(genealogy(:,7)),1);
  end
  totalNumberOfMutations(sample,iteration)=...
    sum(genealogy(:,5))+sum(genealogy(:,6));
end
end
%% Plots and theo

for i=1:size(n,2)
pm=probOfObservingMutation(n(i),100,theta);
subplot(1,2,i)
hold on
 histogram(totalNumberOfMutations)
 plot(0:length(pm)-1,pm*nIterations)
end