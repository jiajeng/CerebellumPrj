function matlabbatch = create_VOI(filepath,matlabbatch,varargin)
in = finputcheck(varargin, ...
       {'TR'        'real'      []  1.5;    
        'roi_name'  'string'    []  []; % roi name 
        'roi_type'  'string'    {'sphere','box'} 'box';
        'cen_cord'  'real'      []  []; % roi center coordinate 
        'dim'       'real'      []  [6 6 6]; % roi size for box
        'rad'       'real'      []  []; % roi size for sphere radius
        'procsNum'  'real'      []  1;
       });

matlabbatch{in.procsNum}.spm.util.voi.spmmat = {fullfile(filepath,'SPM.mat')};
matlabbatch{in.procsNum}.spm.util.voi.adjust = 0;
matlabbatch{in.procsNum}.spm.util.voi.session = 1;
matlabbatch{in.procsNum}.spm.util.voi.name = in.roi_name;
switch in.roi_type
    case 'box'
        matlabbatch{in.procsNum}.spm.util.voi.roi{1}.box.centre = in.cen_cord;
        matlabbatch{in.procsNum}.spm.util.voi.roi{1}.box.dim = in.dim;
        matlabbatch{in.procsNum}.spm.util.voi.roi{1}.box.move.fixed = 1;
    case 'sphere'
        matlabbatch{in.procsNum}.spm.util.voi.roi{1}.sphere.centre = in.cen_cord;
        matlabbatch{in.procsNum}.spm.util.voi.roi{1}.sphere.radius = in.rad;
        matlabbatch{in.procsNum}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
end
matlabbatch{in.procsNum}.spm.util.voi.expression = 'i1';

end




