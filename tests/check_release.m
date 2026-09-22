function check_release()
% Data integrity and small dependency checks, not a full model reproduction.
root=fileparts(fileparts(mfilename('fullpath')));
f=fullfile(root,'02_Data','L600_Data_Consolidated_EN_V02.xlsx');
[w,y,x,ids]=load_l600_submission_ssc(f,'SSC_model_spectra');
assert(isequal(size(x),[372 151]) && isequal(w,(650:2:950)'));
assert(all(isfinite(x),'all') && all(isfinite(y)) && numel(unique(ids))==372);
[wt,g,~]=load_l600_submission_thickness(f,'Thickness_spectra',{'Original peel slice','1 mm peel','5 mm peel'});
assert(isequal(w,wt) && isequal(cellfun(@(a) size(a,1),g),[99 150 142]));
assert(all(cellfun(@(a) all(isfinite(a),'all'),g)));
required={'spxy','projections_qr','airPLS','carspls','randomfrog_pls','fitrsvm','sgolay','pca','TreeBagger'};
for k=1:numel(required)
    assert(exist(required{k},'file')~=0,'Missing dependency: %s. See docs/DEPENDENCIES.md.',required{k});
end
rng(123,'twister');
tr=spxy(x,y,298); assert(numel(unique(tr))==298);
mask=w>=748 & w<=826; a=x(1:30,mask); b=y(1:30);
a=(a-mean(a,2))./std(a,0,2);
p=projections_qr(a,1,5); assert(numel(unique(p))==5);
[c,z]=airPLS(a,1e8,2,0.1,0.05,20); assert(all(isfinite(c),'all') && isequal(size(z),size(a)));
cars=carspls(a,b,3,3,'center',8,0,0,1); assert(~isempty(cars.vsel));
frog=randomfrog_pls(a,b,3,'center',15,2,'regcoef'); assert(isstruct(frog));
m=fitrsvm(a,b,'KernelFunction','gaussian','Standardize',true);
assert(all(isfinite(predict(m,a))));
fprintf('PASS: 391 thickness records, 372 SSC records, 151 wavelengths; dependencies and small-fit checks.\n');
end
