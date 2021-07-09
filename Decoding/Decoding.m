function out=Decoding(decdata,opt)
%%%
% Input:
%     decdata is structure, including two structures (train and test, each including two 4-D matrix s1 and s2)
%     decdata.train.s1 = (cellID * trlN * bins * rpts) matrix 
%%%
arguments
    decdata (1,1) struct
    opt.decoder (1,:) char {mustBeMember(opt.decoder,{'SVM','LDA','NB'})} ='SVM';
end
    
bins=size(decdata.train.s1,3);
trlN=size(decdata.train.s1,2);
rpts=size(decdata.train.s1,4);

cv=cvpartition(trlN,'KFold',10);
y=[zeros(trlN,1);ones(trlN,1)];

out=struct();
cvcorr=cell(0);
shufcorr=cell(0);
%     errcorr{bin+3}=[];
tic
cv_trainingID=cell2mat(arrayfun(@(x)training(cv,x),[1:10],'UniformOutput',false));
cv_testID=cell2mat(arrayfun(@(x)test(cv,x),[1:10],'UniformOutput',false));
% s1kf=arrayfun(@(x)decdata.train.s1(:,x,1,1),cv_trainingID,'UniformOutput',false);
% s1kf_1=reshape(s1kf(cv_trainingID),[],10);

for rpt=1:rpts
%     fprintf('rpts %d of %d\n',rpt,rpts);
    for bin=1:bins
        cvresult=[];
        cvshufresult=[];        
            for kf=1:cv.NumTestSets                
                s1kf=decdata.train.s1(:,cv_trainingID(:,kf),bin,rpt);
                s2kf=decdata.train.s2(:,cv_trainingID(:,kf),bin,rpt);
                varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
                s1kf=s1kf(varsel,:);
                s2kf=s2kf(varsel,:);
                %             Xerr=[s1err(:,varsel);s2err(:,varsel)];
                %             yerr=[zeros(size(s1err,1),1);ones(size(s2err,1),1)];
                Xkf=cat(2,s1kf,s2kf)';
                ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
                
                s1Tkf=decdata.test.s1(varsel,cv_testID(:,kf),bin,rpt);
                s2Tkf=decdata.test.s2(varsel,cv_testID(:,kf),bin,rpt);
                XTkf=cat(2,s1Tkf,s2Tkf)';
                yTkf=y([cv_testID(:,kf);cv_testID(:,kf)]);
                yshufTkf=yTkf(randperm(numel(yTkf)));
                
                
                if strcmp(opt.decoder,'SVM')
                    SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);
                    modelPredict=SVMM.predict(XTkf);
                elseif strcmp(opt.decoder,'LDA')
                    LDAM=fitcdiscr(Xkf,ykf);
                    modelPredict=LDAM.predict(XTkf);
                elseif strcmp(opt.decoder,'NB')
                    NBM=fitcnb(Xkf,ykf);
                    modelPredict=NBM.predict(XTkf);
                end
                
                cvresult=cat(1,cvresult,modelPredict==yTkf);
                cvshufresult=cat(1,cvshufresult,modelPredict==yshufTkf);
                %             cv_err_result=CVLDAModel.predict(Xerr)==yerr;
                
                               
                %             errcorr{bin+3}=[errcorr{bin+3};cv_err_result];
            end
            cvcorr{rpt}{bin}=cvresult;
            shufcorr{rpt}{bin}=cvshufresult; 
    end
end
%
% for rpt=1:rpts
% %     fprintf('rpts %d of %d\n',rpt,rpts);
%     for bin=1:bins
%         cvresult=[];
%         cvshufresult=[];
%             for kf=1:cv.NumTestSets
%                 s1kf=decdata.train.s1(:,cv_trainingID(:,kf),bin,rpt);
%                 s2kf=decdata.train.s2(:,cv_trainingID(:,kf),bin,rpt);
%                 varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
%                 s1kf=s1kf(varsel,:);
%                 s2kf=s2kf(varsel,:);
%                 %             Xerr=[s1err(:,varsel);s2err(:,varsel)];
%                 %             yerr=[zeros(size(s1err,1),1);ones(size(s2err,1),1)];
%                 Xkf=cat(2,s1kf,s2kf)';
%                 ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
%                 
%                 s1Tkf=decdata.test.s1(varsel,cv_testID(:,kf),bin,rpt);
%                 s2Tkf=decdata.test.s2(varsel,cv_testID(:,kf),bin,rpt);
%                 XTkf=cat(2,s1Tkf,s2Tkf)';
%                 yTkf=y([cv_testID(:,kf);cv_testID(:,kf)]);
%                 yshufTkf=yTkf(randperm(numel(yTkf)));
%                 
%                 
%                 if strcmp(opt.decoder,'SVM')
%                     SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);
%                     modelPredict=SVMM.predict(XTkf);
%                 elseif strcmp(opt.decoder,'LDA')
%                     LDAM=fitcdiscr(Xkf,ykf);
%                     modelPredict=LDAM.predict(XTkf);
%                 elseif strcmp(opt.decoder,'NB')
%                     NBM=fitcnb(Xkf,ykf);
%                     modelPredict=NBM.predict(XTkf);
%                 end
%                 
%                 cvresult=cat(1,cvresult,modelPredict==yTkf);
%                 cvshufresult=cat(1,cvshufresult,modelPredict==yshufTkf);
%                 %             cv_err_result=CVLDAModel.predict(Xerr)==yerr;
%                 
%                                
%                 %             errcorr{bin+3}=[errcorr{bin+3};cv_err_result];
%             end
%             cvcorr{rpt}{bin}=cvresult;
%             shufcorr{rpt}{bin}=cvshufresult; 
%     end
% end
toc
out.cvcorr=cell2mat(cellfun(@cell2mat,cvcorr','UniformOutput',false));
out.shufcorr=cell2mat(cellfun(@cell2mat,shufcorr','UniformOutput',false));
end