
clear, clc, clf

%rng default
%Inital values
n=100;            %sample size
theta=5; %mutation rate
N=1000;          %Pop. size
u=theta./(2*N);
I=1000;           %Iteration size
nAlleles=zeros(length(theta),I);

tic
for t=1:length(theta)   %For all different values of theta
    for k=1:I %For all iterations
        nodes=n;                    %Number of nodes in tree
        lO=[1:n]';                  %"List of orphans"
        tP=0;                       %Time into the past, used to calculate time
        
        %Zero vectors to used later
        lambda=zeros(n-1,1);        %Coeff. for exp. dist
        G=zeros(2*n-1,7);           %Genealogy. Result-matrix
        
        %Define T~exp(lambda) for each time step
        for i=2:n
            lambda(i)=nchoosek(i,2)/N;
        end
        T=exprnd(1./(lambda(2:end)));
        T=sort(T);% sort the vector from smallest to largest
        
        %Fill in geneology matrix
        for i=1:n-1 % generating MRCA from n samples from population
            if size(lO,1)>1
                CE=coalescenseOrphans(lO); %Coalescent Event. Retrieve the randomly picked orphans
                
                G(lO(CE(1)),3)=nodes+1;   %In column 3 on the row corresponding to the randomly picked orphan...
                G(lO(CE(2)),3)=nodes+1;   %... num of nodes increases with 1
                G(nodes+1,1)=lO(CE(1));   %Store the "orphans" in G
                G(nodes+1,2)=lO(CE(2));
                
                tP=tP+T(i);               %Time into the past. Total time for tree
                G(nodes+1,4)=tP;          %Store tP in G
                
                nodes=nodes+1;            %Num of nodes increases for every loop
                
                %Update the list of possible descendants such that those used can't be
                %used again
                lO(CE(1))=[];
                if CE(1)<CE(2)
                    lO(CE(2)-1)=[];
                else
                    lO(CE(2))=[];
                end
                lO=[lO;nodes];
            else
                
            end
        end
        
        %------Add mutations to the tree------%
        T=sort(T,'descend');
        
        for i=n+1:2*n-1 %Exclude the first nodes without descendants
            for j=1:2   %For the two different descendants
                timeToDiscard=G(G(i,j),4); %Calculate the time to the most recent descendant
                beta=u(t)*(G(i,4)-timeToDiscard);%Coeff. for Poisson, u times length of branch
                G(i,4+j)=poissrnd(beta); %Take random number from Pois-dist and insert in G
            end
        end
        
        nG=0; %new genes
        for i=0:2*n-2
            for j=1:2 %For the two descendants
                nG=nG+1; %number of different genes increases for each mutation
                if G(2*n-1-i,4+j)>0  %If mutation has occured
                    G(G(2*n-1-i,j),7)=nG; %Add num. of new genes
                elseif G(2*n-1-i,j)==0 %If no mutation occured
                else                   %Do nothing
                    G(G(2*n-1-i,j),7)=G(2*n-1-i,7);
                end
            end
        end
        
        nAlleles(t,k)=size(unique(G(:,7)),1); %Number of different alleles
        nMutations(t,k)=sum(G(:,5))+sum(G(:,6)); %Number of mutations
    end
   
end
toc

%--------- Theory ----------%

j_max=70;
x=probOfObservingMutation(n,j_max,theta);

%--------- Plot the result ----------%
histogram(nMutations,'Facecolor','yellow')
hold on
plot(1:length(x),I.*x,'green','Linewidth',2);
legend('Simulation, n=100','Theory')
