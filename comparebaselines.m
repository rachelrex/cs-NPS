numpulses=length(pulses(:,1));
index=linspace(1,numpulses,numpulses)';

newsizeA=newsize(:,1);
newsizeB=newsize(:,2);
[sortedA,I]=sort(newsizeA);
sortedB=newsizeB(I);
sortednewsize=[sortedA sortedB];

%%
figure
scatter(index,sortedA)
hold on
scatter(index,sortedB)
xlabel('Index')
ylabel('Original Calculated Diameter (microns)')
title('AL 313')
legend('Original Analysis','New Baseline Analysis')

%%
differencevect=sortedB-sortedA;
figure
scatter(sortedA,differencevect)
xlabel('Original Calculated Diameter (microns)')
ylabel('New Diameter - Original Diameter')
title('AL 313')