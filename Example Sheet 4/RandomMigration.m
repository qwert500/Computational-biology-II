function u=RandomMigration(d,m)
u=zeros(m*size(d,1),1);

while size(unique(u),1)~= m*size(d,1)
    
    u=randi([1 size(d,1)], m*size(d,1),1);
    
end
end