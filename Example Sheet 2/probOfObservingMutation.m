function probOfObservingMutation=probOfObservingMutation(sampleSize,maximumNumberOfMutations,theta)
% Task 2 comp bio 2, recursion formula
pPreviuos=zeros(maximumNumberOfMutations+1,1);


for i=2:sampleSize % i=n=sample size
  pCurrent=zeros(maximumNumberOfMutations+1,1);
  for j=0:maximumNumberOfMutations
    if i==2
      pCurrent(j+1)=(theta/(i-1+theta))^j*(i-1)/(i-1+theta);
    else
      for l=0:j
      pCurrent(j+1)=pCurrent(j+1)+pPreviuos(j+1-l)*...
        (theta/(i-1+theta))^l*((i-1)/(i-1+theta));
      end
    end
  end
  pPreviuos=pCurrent;
end
probOfObservingMutation=pCurrent;
end