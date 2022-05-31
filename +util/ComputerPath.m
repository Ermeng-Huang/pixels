function [CodePath,HomePath]=ComputerPath(task)
if ispc
    switch (getenv('USERNAME'))
        case 'hem'
            if strcmp(task,'AIopto')
            CodePath='D:\code\AIopto';
            HomePath='I:\pixel-optogenetic';
            else
                
            end
        case 'zx'
            CodePath='D:\code-hem';
            if strcmp(task,'AIopto')
                HomePath='D:\pixel-optogenetic';
            elseif strcmp(task,'dualtask')
                HomePath='F:\pixel-dualtask';
            end
    end
else
    CodePath='/home/hem/Code';
    HomePath='/media/HDD0/hem/datashare/AI-opto';
end
end