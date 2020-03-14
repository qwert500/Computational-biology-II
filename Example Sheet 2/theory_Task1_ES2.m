%--------Theory--------%
load stirlingsNumber.mat %Taken from Mathematica
S=stirlingsNumber;
m=1:n; %number of possible allelic types

%Ewens sampling formula
nominator=zeros(length(theta),n);
denominator=zeros(length(theta),n-1);
P_n=zeros(length(theta),n);
K=zeros(length(theta),1);
for p=1:size(theta,2) %for all different theta
    for k=2:n
        denominator(p,k-1)=(k-1)+theta(p);
    %Calculate theoretical expected value
        K(p)=K(p) + theta(p)./((k-1)+theta(p));
    end

    nominator(p,:)=abs(S).*theta(p).^(m-1);    
    P_n(p,:)=nominator(p,:)./prod(denominator(p,:),2);
end