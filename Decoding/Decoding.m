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


y=[zeros(trlN,1);ones(trlN,1)];
out=struct();
out.cvcorr=cell(1,bins);
out.shufcorr=cell(1,bins);
cv=cvpartition(trlN,'KFold',10);

cv_trainingID=cell2mat(arrayfun(@(x)training(cv,x),[1:10],'UniformOutput',false));
cv_testID=cell2mat(arrayfun(@(x)test(cv,x),[1:10],'UniformOutput',false));
% s1kf=arrayfun(@(x)decdata.train.s1(:,x,1,1),cv_trainingID,'UniformOutput',false);
% s1kf_1=reshape(s1kf(cv_trainingID),[],10);
switch(nnz(contains(fieldnames(decdata),'test')))
    case 1 
        for rpt=1:rpts
            fprintf('rpts %d of %d\n',rpt,rpts);
            for bin=1:bins                
                for kf=1:cv.NumTestSets
                    s1kf=decdata.train.s1(:,cv_trainingID(:,kf),bin,rpt);
                    s2kf=decdata.train.s2(:,cv_trainingID(:,kf),bin,rpt);
                    varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
                    s1kf=s1kf(varsel,:);
                    s2kf=s2kf(varsel,:);                    
                    Xkf=cat(2,s1kf,s2kf)';
                    ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
                    
                    s1Tkf=decdata.test.s1(varsel,cv_testID(:,kf),bin,rpt);
                    s2Tkf=decdata.test.s2(varsel,cv_testID(:,kf),bin,rpt);
                    XTkf=cat(2,s1Tkf,s2Tkf)';
                    yTkf=y([cv_testID(:,kf);cv_testID(:,kf)]);
                    
%                     yshufTkf=yTkf(randperm(numel(yTkf)));
                    %             Xerr=[s1err(:,varsel);s2err(:,varsel)];
                    %             yerr=[zeros(size(s1err,1),1);ones(size(s2err,1),1)];
                    
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
                    
                    out.cvcorr{:,bin}=cat(1,out.cvcorr{:,bin},modelPredict==yTkf);
                     for i=1:50
                        yshufTkf=yTkf(randperm(numel(yTkf)));
                        cvshufresult=modelPredict==yshufTkf;
                        out.shufcorr{:,bin}=cat(1,out.shufcorr{:,bin},cvshufresult);
                    end   
%                     out.shufcorr{:,bin}=cat(1,out.shufcorr{:,bin},modelPredict==yshufTkf);
                    %             cv_err_result=CVLDAModel.predict(Xerr)==yerr;
                    
                    
                    %             errcorr{bin+3}=[errcorr{bin+3};cv_err_result];
                end
                
            end
        end
    case 2   
        out.cvcorr2=cell(1,bins);
        for rpt=1:rpts
            fprintf('rpts %d of %d\n',rpt,rpts);
            for bin=1:bins                
                for kf=1:cv.NumTestSets
                    s1kf=decdata.train.s1(:,cv_trainingID(:,kf),bin,rpt);
                    s2kf=decdata.train.s2(:,cv_trainingID(:,kf),bin,rpt);
                    varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
                    s1kf=s1kf(varsel,:);
                    s2kf=s2kf(varsel,:);                    
                    Xkf=cat(2,s1kf,s2kf)';
                    ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
                    
                    s1Tkf=decdata.test1.s1(varsel,cv_testID(:,kf),bin,rpt);
                    s2Tkf=decdata.test1.s2(varsel,cv_testID(:,kf),bin,rpt);
                    XTkf1=cat(2,s1Tkf,s2Tkf)';
                    yTkf1=y([cv_testID(:,kf);cv_testID(:,kf)]);
                    
                    
                    XTkf2=[decdata.test2.s1(varsel,:,bin,rpt),decdata.test2.s2(varsel,:,bin,rpt)]';
                    yTkf2=[zeros(size(decdata.test2.s1,2),1);ones(size(decdata.test2.s1,2),1)];
                    
                    if strcmp(opt.decoder,'SVM')
                        SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);
                        modelPredict=SVMM.predict(XTkf1);
                    elseif strcmp(opt.decoder,'LDA')
                        LDAM=fitcdiscr(Xkf,ykf);
                        modelPredict=LDAM.predict(XTkf1);
                    elseif strcmp(opt.decoder,'NB')
                        NBM=fitcnb(Xkf,ykf);
                        modelPredict=NBM.predict(XTkf1);
                    end
                    
                    out.cvcorr{:,bin}=cat(1,out.cvcorr{:,bin},modelPredict==yTkf1);
                    out.cvcorr2{:,bin}=cat(1,out.cvcorr2{:,bin},SVMM.predict(XTkf2)==yTkf2);
                    
                    for i=1:200
                        yshufTkf=yTkf1(randperm(numel(yTkf1)));
                        cvshufresult=modelPredict==yshufTkf;
                        out.shufcorr{:,bin}=cat(1,out.shufcorr{:,bin},cvshufresult);
                    end                    
                end               
               
            end
        end
end
for bin=1:bins  
    [p_temp,~,~]=statistics.permutationTest(out.cvcorr{:,bin},out.shufcorr{:,bin},500);    
    if p_temp<0.001/(rpts*bins)
        out.p(:,bin)=1;
    else
        out.p(:,bin)=0;
    end   
end
end

