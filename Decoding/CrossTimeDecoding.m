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

% parpool(4)
parfor rpt=1:rpts
    cvresult=cell(bins,bins);
    cvshufresult=cell(bins,bins);
    cv_trainingID=cell2mat(arrayfun(@(x)training(cv,x),1:10,'UniformOutput',false));
    cv_testID=cell2mat(arrayfun(@(x)test(cv,x),1:10,'UniformOutput',false));
    
    fprintf('rpts %d of %d\n',rpt,rpts);
    for bin=1:bins       
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
            yshufkf=ykf(randperm(numel(ykf)));
            if strcmp(opt.decoder,'SVM')
                SVMM=fitcsvm(Xkf,ykf,'KernelFunction','linear','Standardize',true);         
                SVMMshuf=fitcsvm(Xkf,yshufkf,'KernelFunction','linear','Standardize',true); 
            elseif strcmp(opt.decoder,'LDA')
                LDAM=fitcdiscr(Xkf,ykf);               
            elseif strcmp(opt.decoder,'NB')
                NBM=fitcnb(Xkf,ykf);
            end
            for bin2=1:bins
                s1Tkf=decdata.test.s1(varsel,cv_testID(:,kf),bin2,rpt);
                s2Tkf=decdata.test.s2(varsel,cv_testID(:,kf),bin2,rpt);
                XTkf=cat(2,s1Tkf,s2Tkf)';
                yTkf=y([cv_testID(:,kf);cv_testID(:,kf)]);                
                                           
                cvresult{bin,bin2}=SVMM.predict(XTkf)==yTkf;
                cvshufresult{bin,bin2}=SVMMshuf.predict(XTkf)==yTkf;
                
               
                
%                 for i=1:50
%                     yshufTkf=yTkf(randperm(numel(yTkf)));
%                     cvshufresult=modelPredict==yshufTkf;
%                     shufcorr{rpt}{bin,bin2}=cvshufresult;
%                 end
                
                %             cv_err_result=CVLDAModel.predict(Xerr)==yerr;              
                %             errcorr{bin+3}=[errcorr{bin+3};cv_err_result];                
            end
   
        end
    end
    cvcorr{rpt}=cellfun(@(x)nnz(x)/length(x),cvresult);
    shufcorr{rpt}=cellfun(@(x)nnz(x)/length(x),cvshufresult);
end

for bin=1:bins
     for bin2=1:bins
         out.cvcorr{bin,bin2}=cellfun(@(x)x(bin,bin2),cvcorr);
         out.shufcorr{bin,bin2}=cellfun(@(x)x(bin,bin2),shufcorr);
     end
end

out.P=cellfun(@(x,y)ranksum(x,y),out.cvcorr,out.shufcorr);


out.p_permutation=cellfun(@(x,y)statistics.permutationTest(x,y,100),out.cvcorr,out.shufcorr);
% out.P=arrayfun(@(x)x<0.001/(rpts*bins),p_temp);


end