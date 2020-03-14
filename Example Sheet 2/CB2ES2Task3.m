clear, clc
%--------------------------
%       Computational
%        Biology II
%     Examples sheet 2
%          Task 3 c)
%--------------------------

%--------------------------
%     Initialization
%--------------------------
n=50;%sample size
theta=100;
N=1000;
%tau=[0.1*N N 10*N]; % all cases
%Nb=[0.5*N 0.1*N 0.01*N]; % all cases
tau=10*N;
Nb=0.01*N;
nIterations=1000;

u=theta./(2*N);
nDiffrentAlleles=zeros(size(theta,2),nIterations);
nMutations=zeros(size(theta,2),nIterations);
T=zeros(n-1,1);
nSingletons=zeros(size(0.1:0.1:5,2),nIterations);
sbIndex=0;
tic
for sb=0.1:0.1:5
  sb
  sbIndex=1+sbIndex;
  taub=Nb*sb; % one Nb for every run of the program
  for nTheta=1:size(theta,2)
    for iteration=1:nIterations
      coalescentEvent=n;
      listOfOrphans=[1:n]';% initial list of ancestorless populations
      timeIntoThePast=0;                       
      
      %-----------------------
      %    Reset parameters
      %   for next iteration
      %-----------------------
      genealogy=zeros(2*n-1,7);           
      
      %-----------------------
      %      Generating
      %   Coalescense time
      %-----------------------
      NCurrent=N;
      T=zeros(n-1,1);
      indexTaub=0;
      indexTaubtau=0;
      for j=2:n
        i=n+2-j;
        lambda=nchoosek(i,2)/NCurrent;
        T(j-1)=exprnd(1/lambda);
        if sum(T)>tau && sum(T)<tau+taub && indexTaub==0
          NCurrent=Nb;
          indexTaub=j;
        elseif sum(T)>tau+taub && indexTaub>0 && indexTaubtau==0
          NCurrent=N;
          indexTaubtau=j;
        end
      end
      
      if indexTaub==0 && indexTaubtau==0
        T=sort(T,'ascend');
      elseif indexTaubtau==0 && indexTaub>0
        T(1:indexTaub-1)=sort(T(1:indexTaub-1),'ascend');
        T(indexTaub:end)=sort(T(indexTaub:end),'ascend');
      elseif indexTaub>0 && indexTaubtau>0
        T(1:indexTaub-1)=sort(T(1:indexTaub-1),'ascend');
        T(indexTaub:indexTaubtau-1)=...
          sort(T(indexTaub:indexTaubtau-1),'ascend');
        T(indexTaubtau:end)=sort(T(indexTaubtau:end),'ascend');
      else
        print('ERROR in time sorting')
      end

      %--------------------------
      %     Generate MRCA
      %    and filling in
      %    genealogy table
      %--------------------------
      for i=1:n-1
        if size(listOfOrphans,1)>1
          orphan=coalescenseOrphans(listOfOrphans); 
          
          genealogy(listOfOrphans(orphan(1)),3)=coalescentEvent+1;   
          genealogy(listOfOrphans(orphan(2)),3)=coalescentEvent+1;
          genealogy(coalescentEvent+1,1)=listOfOrphans(orphan(1));   
          genealogy(coalescentEvent+1,2)=listOfOrphans(orphan(2));
          
          timeIntoThePast=timeIntoThePast+T(i);
          genealogy(coalescentEvent+1,4)=timeIntoThePast;
          
          coalescentEvent=coalescentEvent+1; % New node for next ancestor
          
          
          listOfOrphans(orphan(1))=[]; % Deleting former orphans
          if orphan(1)<orphan(2)
            listOfOrphans(orphan(2)-1)=[];
          else
            listOfOrphans(orphan(2))=[];
          end
          listOfOrphans=[listOfOrphans;coalescentEvent]; 
        else
          
        end
      end
      
      %----------------------
      %       Mutation
      %----------------------
      T=fliplr(T);
      for i=n+1:2*n-1 %Exclude the first nodes without descendants
        for j=1:2   %For the two different descendants
          timeToDiscard=genealogy(genealogy(i,j),4); 
          beta=u(nTheta)*(genealogy(i,4)-timeToDiscard);
          genealogy(i,4+j)=poissrnd(beta); 
        end
      end
      %------------------------
      %   Numrating alleles
      %    for population
      %------------------------
      nG=0; %new genes
      for i=0:2*n-2
        for j=1:2 %For the two descendants
          nG=nG+1; %number of different genes increases for each mutation
          if genealogy(2*n-1-i,4+j)>0  %If mutation has occured
            genealogy(genealogy(2*n-1-i,j),7)=nG; %Add num. of new genes
          elseif genealogy(2*n-1-i,j)==0 %If no mutation occured
          else                   %Do nothing
            genealogy(genealogy(2*n-1-i,j),7)=genealogy(2*n-1-i,7);
          end
        end
      end
      %----------------------
      %   Collecting data
      %----------------------
      nDiffrentAlleles(nTheta,iteration)=size(unique(genealogy(:,7)),1); 
      nMutations(nTheta,iteration)=sum(genealogy(:,5))+sum(genealogy(:,6)); 
      nSingletons(sbIndex,iteration)=...
        FindingNumberOfSingletons(genealogy,n);
    end
  end
end
toc
%% Plots

sb=0.1:0.1:5;

plot(sb,mean(nSingletons,2))









