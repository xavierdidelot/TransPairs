clear all
close all
treefile='trees.nwk';
nrep=1;
for rep=1:nrep
    rep
    unix(sprintf('sed -n %dp %s > tmp.nwk',rep,treefile));
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
            res(i,j)=likelihoodSEIR(tij,ti,tj);
            %npoints=10;
            %for si=1:npoints
            %    s=tij+(min(ti,tj)-tij)*si/(npoints+1);
            %    res(i,j)=res(i,j)+exppdf(abs(tij-s),1/neg)*exppdf(abs(tj-s),1/gamma);
            %end
            %res(i,j)=res(i,j)*(1-expcdf(abs(tij-ti),1/gamma));
        end
    end
    allres=allres+res/nrep;
end
%allres=allres/max(max(allres));%make max=1
allres=log10(allres);

subplot(1,2,1);
plot(tree);
subplot(1,2,2);
imagesc(log10(res));
set(gcf,'Color','w');

f=fopen('matrix.csv','w');
for i=1:n
    fprintf(f,'%f',allres(i,1));
    for j=2:n
        fprintf(f,',%f',allres(i,j));
    end
    fprintf(f,'\n');
end
fclose(f);