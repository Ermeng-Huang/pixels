function out=CrossTimeDecoding(decdata,opt)
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
out.cvcorr=cell(bins,bins);
out.shufcorr=cell(bins,bins);
%     errcorr{bin+3}=[];

cv_trainingID=cell2mat(arrayfun(@(x)training(cv,x),[1:10],'UniformOutput',false));
cv_testID=cell2mat(arrayfun(@(x)test(cv,x),[1:10],'UniformOutput',false));

for rpt=1:rpts
    fprintf('rpts %d of %d\n',rpt,rpts);
    for kf=1:cv.NumTestSets
        for bin=1:bins
            for bin2=1:bins            
                s1kf=decdata.train.s1(:,cv_trainingID(:,kf),bin,rpt);
                s2kf=decdata.train.s2(:,cv_trainingID(:,kf),bin,rpt);
                varsel=any(s1kf-s1kf(:,1),2) & any(s2kf-s2kf(:,1),2);
                s1kf=s1kf(varsel,:);
                s2kf=s2kf(varsel,:);
                %             Xerr=[s1err(:,varsel);s2err(:,varsel)];
                %             yerr=[zeros(size(s1err,1),1);ones(size(s2err,1),1)];
                Xkf=cat(2,s1kf,s2kf)';
                ykf=y([cv_trainingID(:,kf);cv_trainingID(:,kf)]);
                
                s1Tkf=decdata.test.s1(varsel,cv_testID(:,kf),bin2,rpt);
                s2Tkf=decdata.test.s2(varsel,cv_testID(:,kf),bin2,rpt);
                XTkf=cat(2,s1Tkf,s2Tkf)';
                yTkf=y([cv_testID(:,kf);cv_testID(:,kf)]);
                
                Train{bin,bin2}=[Xkf,ykf];
                Test{bin,bin2}=[XTkf,yTkf];
            end
        end  
        
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
                
                cvresult=modelPredict==yTkf;
                out.cvcorr{bin,bin2}=cat(1,out.cvcorr{bin,bin2},cvresult);
                
                for i=1:1000
                    yshufTkf=yTkf(randperm(numel(yTkf)));
                    cvshufresult=modelPredict==yshufTkf;
                    out.shufcorr{bin,bin2}=cat(1,out.shufcorr{bin,bin2},cvshufresult);
                end
                
                %             cv_err_result=CVLDAModel.predict(Xerr)==yerr;              
                %             errcorr{bin+3}=[errcorr{bin+3};cv_err_result];
                

    end
end

for bin=1:bins
    for bin2=1:bins
        [p_temp,~,~]=statistics.permutation(out.corr,out.shufcorr,1000);
        [p_temp,~,~]=signrank(out.cvcorr{bin,bin2},out.shufcorr{bin,bin2});
        if p_temp<0.001/(rpts*bins)
            out.p(bin,bin2)=1;
        else
            out.p(bin,bin2)=0;
        end
    end
end

end