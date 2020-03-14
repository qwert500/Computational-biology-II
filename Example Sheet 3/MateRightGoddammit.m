function mateRight=MateRightGoddammit(diploid,n)

condition=0;
while condition==0
    randomNumber=coalescenseOrphans(diploid);
    newDiploidCoalesce=diploid(randomNumber(1),:)+diploid(randomNumber(2),:);
    
    %if newDiploidCoalesce(:,1)==0 || newDiploidCoalesce(:,2)==0
    %elseif newDiploidCoalesce(:,1)>n || newDiploidCoalesce(:,2)>n
    %else
        condition=1;
    %end
end

mateRight=randomNumber;


end