%-------------------------
%      Computational
%       Biology II
%         Task 3
%-------------------------
clc, clear
%-------------------------
%     Initialization
%-------------------------
ns=[2 4 8 16];
Rs=[0 10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 1 2 3 5 10 20];
nIterations=3000;

nGenes=2;
tau=zeros(2,nIterations,size(Rs,2),size(ns,2));
timeIntoThePast=0;
%-------------------------
%     algorithm
%-------------------------
tic
for sample=1:1:size(ns,2)
  n=ns(sample);
  for nR=1:size(Rs,2)
    R=Rs(nR);
    for iteration=1:1:nIterations
      numberOfLinages=n;
      timeIntoThePast=0;
      diploid=ones(n,2);
      while tau(1, iteration,nR,sample)==0 || tau(2,iteration,nR,sample)==0
        
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
        
        %-------------------------
        %     Time to MRCA
        %-------------------------
        
        for i=1:size(diploid,1) % Check if MRCA is found for any gene and sets the time
          for j=1:nGenes
            if diploid(i,j)==n && tau(j,iteration,nR,sample)==0
              
              tau(j,iteration,nR,sample)=timeIntoThePast;
            end
          end
        end
        
        numberOfLinages=size(diploid,1);
      end
    end
  end
end
toc
%%
%-------------------------
%      Computational
%       Biology II
%         Task 3
%-------------------------

%-------------------------
%    Initialization
%-------------------------
tauA=zeros(nIterations,size(Rs,2),size(ns,2));
tauB=zeros(nIterations,size(Rs,2),size(ns,2));
for sample=1:size(ns,2)
  for nR=1:1:size(Rs,2)
    tauA(:,nR,sample)=tau(1,:,nR,sample);
    tauB(:,nR,sample)=tau(2,:,nR,sample);
    TAU(nR,sample)=mean(tauA(:,nR,sample))./2+mean(tauB(:,nR,sample))./2;
    TAUsquare(nR,sample)=mean(tauA(:,nR,sample).^2)./2+mean(tauB(:,nR,sample).^2)./2;
  end
end
for sample=1:size(ns,2)
  for nR=1:1:size(Rs,2)
    correlationFunction(nR,sample)=(mean(tauA(:,nR,sample).*tauB(:,nR,sample))-mean(tauA(:,nR,sample))*mean(tauB(:,nR,sample)))/...
      (TAUsquare(nR,sample)-TAU(nR,sample)^2);
    theoreticalCorrelationFunction(nR)=(Rs(nR)+18)/(Rs(nR)^2+13*Rs(nR)+18);
  end
end
for sample=2:size(ns,2)
  subplot(1,size(ns,2)-1,sample-1)
  plot(log(Rs),correlationFunction(:,1),'bo--')
  hold on
  plot(log(Rs),theoreticalCorrelationFunction,'g')
  plot(log(Rs),correlationFunction(:,sample),'ro--')
  xlabel('log(R)')
  ylabel('Correlation function')
  title(['Correlation function, n=2, n=',num2str(ns(sample))])
  legend('Simulation n=2','Theory',['simulation n=',num2str(ns(sample))])
end




