function synergy2(datain, dataout)

%%

load('result/expedstringmat.txt')
load(datain)
load('result/expdata.txt')

norin = expedstringmat;
clear expedstringmat;
F0 = genelist;
geneexp = expdata;

inter = 200;
changethreshold = 0.0002;
alpha = 0.7;
 

todelete = find(sum(norin,2) == 0);
norin(todelete,:) = [];
norin(:,todelete) = [];
geneexp(todelete,:) = [];
F0(todelete,:) = [];

norA = sqrt(1 ./ sum(norin,2));
norout = norA.*norin.*norA';
norout = geneexp.*norout.*geneexp';
A=norout;
clear norout;

Ft = F0;
for n = 1:inter
	Ft2 = alpha*A*Ft + (1-alpha)*F0;
	change = sum(abs(Ft2 - Ft));
    if(change < changethreshold)
    	break;
    end;
end
writematrix(Ft2,dataout)

