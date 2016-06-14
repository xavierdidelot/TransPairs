clear all
close all

treefile='trees.nwk';%This file contains a sample of trees, for example the output from running BEAST
nrep=10;%Number of trees to be used from the sample

for rep=1:nrep
    rep
    unix(sprintf('sed -n %dp %s > tmp.nwk',rep,treefile));%Extract corresponding line
    tree=phytreeread('tmp.nwk');
    n=get(tree,'NUMLEAVES');
    res=zeros(n,n);
    if rep==1,allres=res;end
    for i=1:n
        for j=1:n
            if i==j,continue;end
            %Consider transmission from i to j
            tree2=prune(tree,setdiff(1:n,[find(getbyname(tree,num2str(i),'EXACT','true')) find(getbyname(tree,num2str(j),'EXACT','true'))]));
            if get(tree2,'NUMLEAVES')~=2,'error',end
            d=get(tree2,'DISTANCES');
            ti=d(find(getbyname(tree2,num2str(i),'EXACT','true')));
            tj=d(find(getbyname(tree2,num2str(j),'EXACT','true')));
            tij=0;
            if (true) 
	        res(i,j)=likelihoodSEIR(tij,ti,tj);%Calculate pairwise likelihood using the likelihoodSEIR function
	    else
	        %This is a previous attempt at calculating the pairwise likelihood by averaging the integral over several points
                npoints=10;
                for si=1:npoints
                    s=tij+(min(ti,tj)-tij)*si/(npoints+1);
                    res(i,j)=res(i,j)+exppdf(abs(tij-s),1/neg)*exppdf(abs(tj-s),1/gamma);
                end
                res(i,j)=res(i,j)*(1-expcdf(abs(tij-ti),1/gamma));
            end
        end
    end
    allres=allres+res/nrep;
end

allres=log10(allres);

%Plot the pairwise likelihood values side-by-side with the phylogeny
if (false) 
    subplot(1,2,1);
    plot(tree);
    subplot(1,2,2);
    imagesc(allres);
    set(gcf,'Color','w');
end

%Output the pairwise likelihood values into the 'matrix.csv' file
f=fopen('matrix.csv','w');
for i=1:n
    fprintf(f,'%f',allres(i,1));
    for j=2:n
        fprintf(f,',%f',allres(i,j));
    end
    fprintf(f,'\n');
end
fclose(f);

%Call R script 'edmunds.R' to find the optimal transmission tree given the pairwise likelihood values
system('Rscript edmunds.R')
