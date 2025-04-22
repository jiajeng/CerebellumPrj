function checkfile(file,filename,filepath)
    if size(file,1)>1
        msg = 'find more than one file: ';
        for i = 1:size(file,1)-1
            msg = [msg,file(i),', '];
        end
        msg = [msg,file(end)];
        error(msg); 
    end
    if isempty(file), error(['can not find ',filename,' in ',filepath]); end
end     
