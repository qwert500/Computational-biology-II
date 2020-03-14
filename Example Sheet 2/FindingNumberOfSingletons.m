function findNumberOfSingletons=FindNumberOfSingletons(genealogy,n)
%alleles is a column
nSingletons=0;
for i=1:n
  for j=n+1:n+n-1
    for individual=1:2
      if i==genealogy(j,individual)
        nSingletons=nSingletons+genealogy(j,4+individual);
      end
    end
  end
end
findNumberOfSingletons=nSingletons;
end