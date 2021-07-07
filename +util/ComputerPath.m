function [CodePath,HomePath]=ComputerPath()
if ispc
    switch (getenv('USERNAME'))
        case 'hem'
            CodePath='D:\code\AIopto';
            HomePath='I:\pixel-optogenetic';
        case 'zx'
            CodePath='D:\code-hem';
            HomePath='D:\pixel-optogenetic';
    end
else
    CodePath='/home/hem/Code';
    HomePath='/media/HDD0/hem/datashare/AI-opto';
end
end