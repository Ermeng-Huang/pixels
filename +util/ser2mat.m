
%% Xiaoxing training box
function mat=ser2mat(file)
% SER2MAT  Convert java recorded serial messages to Matlab matrix.
%   mat = SER2MAT('file') convert a file.
%   @ZX, 2014.10

fis=javaObject('java.io.FileInputStream',file);
ois=javaObject('java.io.ObjectInputStream',fis);
evts=cell2mat(cell(ois.readObject().toArray()));
mat=reshape(evts,3,size(evts,1)/3);


