function findingNumberOfUniqueAlleles=FindingNumberOfUniqueAlleles(allele)
%alleles is a column, unique in sence of that it only exist one in the
%sample
singleton=0;
numberOfAlleles=size(allele,1);
for i=1:numberOfAlleles
if size(find(allele(i)==allele(:)),1)==1
  singleton=singleton+1;
end
end
findingNumberOfUniqueAlleles=singleton;