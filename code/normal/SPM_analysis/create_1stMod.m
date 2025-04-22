function [matlabbatch,procsNum] = create_1stMod(filepath,outputpath,matlabbatch,varargin)
in = finputcheck(varargin, ...
       {'TR'        'real'      []  1.5;    
        'T'         'real'      []  16;
        'T0'        'real'      []  8;
        'units'     'string'    {'scans','secs'} 'secs';
        'filename'  'string'    []  'ds*.nii';
        'condition' 'cell'      []  {'HDHF','HDLF','LDHF','LDLF','ERROR'};
        'onsetname' 'string'    {'SoundOnset','WordOnset'}  'SoundOnset';
        'dur'       'real'      []  0;
        'timModul'  'real'      []  0;
        'condPar'   'cell'      []  {'ResponseTime'};
        'procsNum'  'real'      []  [];
        'analysis'  'string'    []  ' ';
        'voiName'   'string'    []  [];
       });
procsNum = in.procsNum;
% outputpath = fullfile(filepath,'parmModu');
if ~exist(outputpath,'dir'), mkdir(outputpath); end
namePat = split(in.filename,'*');
alldataname = string({dir(filepath).name});
dataname = [];
namePatmL = max(cellfun(@length,namePat));
% ls
for j = 1:length(alldataname)
    f = false;
    dn = char(alldataname(j));
    if length(dn) < namePatmL, continue; end 
    if dn(1:length(namePat{1})) == string(namePat{1}) || dn(end-length(namePat{end}):end) == string(namePat{end})
        f = true;
    end
    dn = dn(length(namePat{1}):end-length(namePat{end}));
    for i = 2:length(namePat)-2
        if dn(1:length(namePat{i})) == namePat{i}, f = f && true; end
        dn = dn(length(namePat{i}:end-length(namePat{end})));
    end
    if f, dataname = cat(1,dataname,char(alldataname(j)));end
end
if size(dataname,1) > 1
    warning('select more than one file.');
    keyboard;
end
scansNum = size(niftiread(fullfile(filepath,dataname))); scansNum = scansNum(end);
Ncondition = length(in.condition);
Npar = length(in.condPar);

matlabbatch{procsNum}.spm.stats.fmri_spec.dir = {outputpath};
matlabbatch{procsNum}.spm.stats.fmri_spec.timing.units = in.units;
matlabbatch{procsNum}.spm.stats.fmri_spec.timing.RT = in.TR;
matlabbatch{procsNum}.spm.stats.fmri_spec.timing.fmri_t = in.T;
matlabbatch{procsNum}.spm.stats.fmri_spec.timing.fmri_t0 = in.T0;

%% define scans file
scansfile = append(repmat({fullfile(filepath,dataname)},scansNum,1),repmat({','},scansNum,1),cellfun(@num2str,num2cell(1:scansNum),'UniformOutput',false)');

matlabbatch{procsNum}.spm.stats.fmri_spec.sess.scans = scansfile;

%% define conditions 
if in.analysis == "PPI" || in.analysis == "ppi"
    ppifilename = {dir(outputpath).name};
    ppifilename = ppifilename(contains(ppifilename,"PPI"));
    ppifilename = ppifilename(contains(ppifilename,[in.voiName,'.mat']));

    rNum = 1;
    for ncond = 1:length(ppifilename)
        load(fullfile(outputpath,char(ppifilename(ncond))))

        % variable PPI .ppi --> interaction 
        %              .Y --> VOI_BOLD 
        %              .P --> x_psych
        reg_name = split(PPI.name,'x');
        matlabbatch{procsNum}.spm.stats.fmri_spec.sess.regress(rNum).name = [reg_name{1},'_inter'];
        matlabbatch{procsNum}.spm.stats.fmri_spec.sess.regress(rNum).val = PPI.ppi; % interaction

        matlabbatch{procsNum}.spm.stats.fmri_spec.sess.regress(rNum+1).name = reg_name{1};
        matlabbatch{procsNum}.spm.stats.fmri_spec.sess.regress(rNum+1).val = PPI.P; % x_
        rNum = rNum+2;
    end

    matlabbatch{procsNum}.spm.stats.fmri_spec.sess.regress(end+1).name= reg_name{2};
    matlabbatch{procsNum}.spm.stats.fmri_spec.sess.regress(end).val = PPI.Y;
else
    load(fullfile(filepath,'..','..','BEHAV','behD.mat'));
    for ncond = 1:Ncondition
        for npar = 1:Npar
            %%
            matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).name = in.condition{ncond};
            %%
            matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).onset = cell2mat([behD.(in.condition{ncond}).(in.onsetname)]);
            %%
            matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).duration = in.dur;
            matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).tmod = in.timModul;
            %%
            matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).orth = 1;
    
            param = cell2mat([behD.(in.condition{ncond}).(in.condPar{npar})]);
            if in.condition{ncond} == "ERROR"
                matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).pmod = struct('name', {}, 'param', {}, 'poly', {});
            else
                parName = char(in.condPar{npar});
                if length(parName) > 3
                    parName = parName(isstrprop(parName,'upper'));
                end
                matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).pmod(npar).name = parName;
                matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).pmod(npar).param = param;
                matlabbatch{procsNum}.spm.stats.fmri_spec.sess.cond(ncond).pmod(npar).poly = 1;
            end
        end
    end
end

matlabbatch{procsNum}.spm.stats.fmri_spec.sess.multi = {''};
% matlabbatch{procsNum}.spm.stats.fmri_spec.sess.multi_reg = {'C:\Users\user\Desktop\jeng\CerebellumPrj\data\SUB048\niifile\WORD\test.txt'};
matlabbatch{procsNum}.spm.stats.fmri_spec.sess.multi_reg = {''};
matlabbatch{procsNum}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{procsNum}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{procsNum}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{procsNum}.spm.stats.fmri_spec.volt = 1;
matlabbatch{procsNum}.spm.stats.fmri_spec.global = 'None';
matlabbatch{procsNum}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{procsNum}.spm.stats.fmri_spec.mask = {''};
matlabbatch{procsNum}.spm.stats.fmri_spec.cvi = 'AR(1)';

procsNum = procsNum+1;
matlabbatch{procsNum}.spm.stats.fmri_est.spmmat = {fullfile(outputpath,'SPM.mat')};
matlabbatch{procsNum}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{procsNum}.spm.stats.fmri_est.method.Classical = 1;
end