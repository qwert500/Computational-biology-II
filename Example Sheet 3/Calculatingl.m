function l=Calculatingl(diploid,n)
subtractionTerm=0;
for i=1:size(diploid,1)
  condition=0;
  for j=1:size(diploid,2)
    if (diploid(i,j)==0 || diploid(i,j)==n) && condition==0
      subtractionTerm=subtractionTerm+1;
      condition=1;
    end
  end
end
l=size(diploid,1)-subtractionTerm;
end