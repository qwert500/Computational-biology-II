function coalescenseOrphans=coalescenseOrphans(listOfOrphans)
orphan1=0;
orphan2=0;
while orphan1==orphan2
orphan1=randi(size(listOfOrphans,1));
orphan2=randi(size(listOfOrphans,1));
if orphan1~=orphan2
  coalescenseOrphans=[orphan1 orphan2];
end

end

end
