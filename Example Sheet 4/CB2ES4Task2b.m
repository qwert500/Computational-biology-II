%-------------------------
%      Comp. Bio II
%     Example sheet 4
%         Task 2b       
%-------------------------
clear, clc
%-------------------------
%     Initialization
%-------------------------
m=0.1;
N=200;
s=0.05; %[0.05 0.1];
timeLimit=500;
nIterations=1000; 
nTimeSteps=3000;

probabilityOfFixationA=0;
probabilityOfFixationa=0;
storex1=zeros(nIterations,nTimeSteps+1); %All freq. indep. of situation
storex2=zeros(nIterations,nTimeSteps+1);


w1AA=1-s; w1Aa=1-s/2; w1aa=1;
w2AA=1;   w2Aa=1-s/2; w2aa=1-s;

j=0;
tic
for i=1:nIterations
    %-------------------------
    %     Creating pop.s
    %-------------------------
    diploid1=zeros(N,2); 
    diploid2=zeros(N,2); 
    
    diploid2(1,1)=1;    
    
    x1=0; x2=1/(2*N);
    
    for t=2:nTimeSteps+1
        %-------------------------
        %     Migration
        %-------------------------
        
        u1=RandomMigration(diploid1,m);
        u2=RandomMigration(diploid2,m);
        
        for migrate=1:size(u1,1)
            migratingDiploids1=diploid1(u1(migrate),:);
            miggratingDiploids2=diploid2(u2(migrate),:);
            
            diploid1(u1(migrate),:)=miggratingDiploids2;
            diploid2(u2(migrate),:)=migratingDiploids1;
        end
        
        %Update all proportions and frequencies, pop1
        fAA_1=size(find(sum(diploid1,2)==2),1)/(2*N); %AA
        fAa_1=size(find(sum(diploid1,2)==1),1)/(2*N); %Aa
        faa_1=size(find(sum(diploid1,2)==0),1)/(2*N); %aa
        
        w_bar_1=fAA_1*w1AA+fAa_1*w1Aa+faa_1*w1aa;
        
        proportionOfAA_1=fAA_1*w1AA/w_bar_1;
        proportionOfAa_1=fAa_1*w1Aa/w_bar_1;
        proportionOfaa_1=faa_1*w1aa/w_bar_1;
        ptot_1=proportionOfAA_1+proportionOfAa_1+proportionOfaa_1;
        
        %Update all proportions and frequencies, pop2
        fAA_2=size(find(sum(diploid2,2)==2),1)/(2*N); %AA
        fAa_2=size(find(sum(diploid2,2)==1),1)/(2*N); %Aa
        faa_2=size(find(sum(diploid2,2)==0),1)/(2*N); %aa
        
        w_bar_2=fAA_2*w2AA+fAa_2*w2Aa+faa_2*w2aa;
        
        proportionOfAA_2=fAA_2*w2AA/w_bar_2;
        proportionOfAa_2=fAa_2*w2Aa/w_bar_2;
        proportionOaa_2=faa_2*w2aa/w_bar_2;
        ptot_2=proportionOfAA_2+proportionOfAa_2+paa_2;
        
        
        u=rand(2*N,1);
        
        %Create new diploids (children)
        for k=1:length(u)
            if (u(k)<=proportionOfAA_1)
                newDiploids(k)=1;
            elseif (proportionOfAA_1<u(k)) && (u(k)<=(proportionOfAa_1+proportionOfAA_1))
                newDiploids(k)=randi([0 1]);
            else
                newDiploids(k)=0;
            end
        end
        
        %Children becomes parents
        diploid1=[newDiploids(1:N)' newDiploids(N+1:end)'];
        
        u=rand(2*N,1);
        
        %Create new diploids (children)
        for k=1:length(u)
            if u(k)<proportionOfAA_2
                newDiploids(k)=1;
            elseif (proportionOfAA_2<=u(k)) && (u(k)<=(proportionOfAa_2+proportionOfAA_2))
                newDiploids(k)=randi([0 1]);
            else
                newDiploids(k)=0;
            end
        end
        
        %Children becomes parents
        diploid2=[newDiploids(1:N)' newDiploids(N+1:end)'];
        
        
        %Calculate freq.
        x1(t)=size(find(diploid1==1),1)/(2*size(diploid1,1));
        x2(t)=size(find(diploid2==1),1)/(2*size(diploid2,1));
        
    end
    %Calculate freq. for all cases
    storex1(i,1:size(x1,2))=x1;
    storex2(i,1:size(x1,2))=x2;
    
    %Calculate freq. for only quasi-ss event
    if storex1(i,timeLimit)~=0
        j=j+1;
        probabilityOfFixationA=probabilityOfFixationA+1/nIterations;
        quasix1(j,1:size(x1,2))=x1;
        quasix2(j,1:size(x1,2))=x2;
    else
        probabilityOfFixationa=probabilityOfFixationa+1/nIterations;
    end
    
end
toc

load gong.mat
sound(y,Fs)

%%
clc, clf

%Summarize result
storex1Average=mean(storex1,1);
storex2Average=mean(storex2,1);

fixationsx1=mean(storex1Average(300:end));
fixationAllx2=mean(storex2Average(300:end));

quasix1Average=mean(quasix1,1);
quasix2Average=mean(quasix2,1);

fixationQuasix1=mean(quasix1Average(300:end));
fixationQuasix2=mean(quasix2Average(300:end));

%Theory
load theory_x1m001s01
load theory_x2m001s01
AAA=theory_x1m001s01;
BBB=theory_x2m001s01

x1_theory=ones(1,length(storex1Average))*AAA;
x2_theory=ones(1,length(storex1Average))*BBB;

%Plot result
t=1:length(storex1Average);

subplot(1,2,1)
plot(t,storex1Average,t,storex2Average)
legend('x1','x2')
title('All events')
grid on

subplot(1,2,2)
plot(t,quasix1Average,t,quasix2Average)
hold on
plot(t,x1_theory,'black--',t,x2_theory,'black--')
legend('x1, simulation','x2, simulation','x1, theory','x2, theory')
title('All quasi-events m=0.01,s=0.1')
grid on
%Compare with theory
diff1=fixationQuasix1/AAA
diff2=fixationQuasix2/BBB
%% Testing
t=1:length(storex1Average);
plot(t,quasix1(1:5,:))

