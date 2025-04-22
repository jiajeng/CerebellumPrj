% change conn project path and name 
connP_filename = 'conn_cerebellum_WORD';
connP_path = fullfile('C:\Users\user\Documents\GitHub\Cerebellumprj\code',connP_filename);

analysis_name = 'gPPI_01';

% check roi name is right
roi_name_path = fullfile(connP_path,connP_filename,'results','firstlevel',analysis_name);
load(fullfile(roi_name_path,'_list_sources.mat')); % variable name = sourcenames
com_sourcenames = sourcenames;
keyboard;
effect_names = {'HDHF','HDLF','LDHF','LDLF'};
con = {[-1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1]};
%%
for nroi = 1:length(com_sourcenames)
    for neft = 1:length(con)
    conn_batch( 'filename',fullfile(connP_path,[connP_filename,'.mat']), ...
        'Results.analysis_number',analysis_name, ...
        'Results.between_subjects.effect_names',   {'AllSubjects'}, ...
        'Results.between_subjects.contrast',[1], ...
        'Results.between_conditions.effect_names', effect_names, ...
        'Results.between_conditions.contrast',con{neft}, ...
        'Results.between_sources.effect_names',com_sourcenames(nroi), ...
        'Results.between_sources.contrast',[1], ...
        'Results.display',false)
    end
end
