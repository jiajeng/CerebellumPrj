function matlabbatch = Estimate(filepath,matlabbatch,varargin)
in = finputcheck(varargin, ...
       {'procsNum'  'real'      []  [];
       });


matlabbatch{in.procsNum}.spm.stats.fmri_est.spmmat = {fullfile(filepath,'SPM.mat')};
matlabbatch{in.procsNum}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{in.procsNum}.spm.stats.fmri_est.method.Classical = 1;


end