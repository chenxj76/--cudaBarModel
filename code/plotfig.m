T=readtable('U300');
area=table2array(T);
imagesc(1:300,1:300,area(1:300,1:300))%eval: Execute MATLAB expression in text
axis([1 300 1 300])
colorbar();