% input data from
clear
opt.pair='ext';
opt.type='dualtask';
xcorr_folder='bzdata_10ms';
homedir='F:\pixel-dualtask';
[sig,pair]=load_sig_pair(fullfile(homedir,'xcorr','bzdata_10ms','bzdata_zx'),'ext');
allneuron=false;
n=4;

sessionIdx=double(unique(floor(sig.suid/100000)));
NeuronSet_sel=util.ChooseNeuronSet_dualtask('sel','sel');
NeuronSet_all=util.ChooseNeuronSet_dualtask('sel','all');
NeuronSet_d2=util.ChooseNeuronSet_dualtask('sel','d2');
NeuronSet_d3=util.ChooseNeuronSet_dualtask('sel','d3');

load(fullfile(homedir,'difftype_ByPEV_1019.mat'))

load(fullfile(homedir,'reactived-neuron.mat'),'NeuronSet_gain');


%% chain_7
chain_num=5;
chain=[];
id_begin=NeuronSet_sel.id;
id_end=NeuronSet_d2.id;
for i=1:length(sessionIdx)
   sig_begin=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i)&ismember(sig.suid(:,1),id_begin),:));
   sig_sess=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i)&ismember(sig.suid(:,1),NeuronSet_sel.id),:)); 
%    sig_sess=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i),:)); 
   for ii=1:size(sig_begin,1)
       id1=sig_begin(ii,1);
       id2=sig_begin(ii,2);
       id3_temp=sig_sess(sig_sess(:,1)==id2,2);
       if isempty(id3_temp)
           continue
       else
           for iii=1:length(id3_temp)
               id3=id3_temp(iii);
               id4_temp=sig_sess(sig_sess(:,1)==id3,2);
               if isempty(id4_temp)
                   continue
               else
                   for iiii=1:length(id4_temp)
                       id4=id4_temp(iiii);
                       id5_temp=sig_sess(sig_sess(:,1)==id4,2);
                       if isempty(id5_temp)
                           continue
                       else
                           if chain_num>5
                               for iiiii=1:length(id5_temp)
                                   id5=id5_temp(iiiii);
                                   id6_temp=sig_sess(sig_sess(:,1)==id5,2);
                                   if isempty(id6_temp)
                                       continue
                                   else
                                       if chain_num>6
                                           for iiiiii=1:length(id6_temp)
                                               id6=id6_temp(iiiiii);
                                               id7_temp=sig_sess(sig_sess(:,1)==id6,2);
                                               if isempty(id7_temp)
                                                   continue
                                               else
                                                   if chain_num>7
                                                       for iiiiiii=1:length(id7_temp)
                                                           id7=id7_temp(iiiiiii);
                                                           id8_temp=sig_sess(sig_sess(:,1)==id7&ismember(sig_sess(:,2),NeuronSet_gain.id),2);
                                                           if isempty(id8_temp)
                                                               continue
                                                           else
                                                               chain(end+1:end+length(id8_temp),:)=cat(2,repmat([id1,id2,id3,id4,id5,id6,id7],length(id8_temp),1),id8_temp);
                                                           end
                                                       end
                                                   else
                                                       % chain_6
                                                       id7_temp=sig_sess(sig_sess(:,1)==id6&ismember(sig_sess(:,2),NeuronSet_gain.id),2);
                                                       if isempty(id7_temp)
                                                           continue
                                                       else
                                                           chain(end+1:end+length(id7_temp),:)=cat(2,repmat([id1,id2,id3,id4,id5,id6],length(id7_temp),1),id7_temp);
                                                       end
                                                   end
                                               end
                                           end
                                       else
                                           % chain_6
                                           id6_temp=sig_sess(sig_sess(:,1)==id5&ismember(sig_sess(:,2),NeuronSet_gain.id),2);
                                           if isempty(id6_temp)
                                               continue
                                           else
                                               chain(end+1:end+length(id6_temp),:)=cat(2,repmat([id1,id2,id3,id4,id5],length(id6_temp),1),id6_temp);
                                           end
                                       end
                                   end
                               end
                           else
                               % chain_5
                               id5_temp=sig_sess(sig_sess(:,1)==id4&ismember(sig_sess(:,2),id_end),2);
                               if isempty(id5_temp)
                                   continue
                               else
                                   chain(end+1:end+length(id5_temp),:)=cat(2,repmat([id1,id2,id3,id4],length(id5_temp),1),id5_temp);
                               end
                           end
                       end
                   end
               end
           end
       end
   end
end
save(fullfile(homedir,'xcorr','chain_5_mem.mat'),'chain','-v7.3')
return
%% chain_5
chain_5=[];
for i=1:length(sessionIdx)
%    sig_begin=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i)&ismember(sig.suid(:,1),NeuronSet_sel.id),:));
   sig_begin=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i)&ismember(sig.suid(:,1),gain_all.id),:));
   sig_sess=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i)&ismember(sig.suid(:,1),NeuronSet_sel.id),:)); 
%    sig_sess=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i),:)); 
   for ii=1:size(sig_begin,1)
       id1=sig_begin(ii,1);
       id2=sig_begin(ii,2);
       id3_temp=sig_sess(sig_sess(:,1)==id2,2);
       if isempty(id3_temp)
           continue
       else
           for iii=1:length(id3_temp)
               id3=id3_temp(iii);
               id4_temp=sig_sess(sig_sess(:,1)==id3,2);
               if isempty(id4_temp)
                   continue
               else
                   for iiii=1:length(id4_temp)
                       id4=id4_temp(iiii);
%                        id5_temp=sig_sess(sig_sess(:,1)==id4&ismember(sig_sess(:,2),setdiff(NeuronSet_all.id,NeuronSet_gain.id)),2);
                       id5_temp=sig_sess(sig_sess(:,1)==id4&ismember(sig_sess(:,2),NeuronSet_gain.id),2);
                       if isempty(id5_temp)
                           continue
                       else
                           chain_5(end+1:end+length(id5_temp),:)=cat(2,repmat([id1,id2,id3,id4],length(id5_temp),1),id5_temp);
                       end
                   end
               end
           end
       end
       
   end
end
save(fullfile(homedir,'xcorr','chain_5_gain(1)-mem(2-4)-reac(5).mat'),'chain_5','-v7.3')

%% chain_6
chain_6=[];
for i=1:length(sessionIdx)
   sig_begin=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i)&ismember(sig.suid(:,1),gain_all.id),:));
   sig_sess=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i)&ismember(sig.suid(:,1),NeuronSet_sel.id),:)); 
%    sig_sess=double(sig.suid(floor(sig.suid(:,1)/100000)==sessionIdx(i),:)); 
   for ii=1:size(sig_begin,1)
       id1=sig_begin(ii,1);
       id2=sig_begin(ii,2);
       id3_temp=sig_sess(sig_sess(:,1)==id2,2);
       if isempty(id3_temp)
           continue
       else
           for iii=1:length(id3_temp)
               id3=id3_temp(iii);
               id4_temp=sig_sess(sig_sess(:,1)==id3,2);
               if isempty(id4_temp)
                   continue
               else
                   for iiii=1:length(id4_temp)
                       id4=id4_temp(iiii);
                       id5_temp=sig_sess(sig_sess(:,1)==id4,2);
                       if isempty(id5_temp)
                           continue
                       else
                            for iiiii=1:length(id5_temp)
                                id5=id5_temp(iiiii);
                                id6_temp=sig_sess(sig_sess(:,1)==id5&ismember(sig_sess(:,2),NeuronSet_gain.id),2);
                                if isempty(id6_temp)
                                    continue
                                else
                                    chain_6(end+1:end+length(id6_temp),:)=cat(2,repmat([id1,id2,id3,id4,id5],length(id6_temp),1),id6_temp);
                                end
                            end
                       end
                   end
               end
           end
       end
       
   end
end
save(fullfile(homedir,'xcorr','chain_6_gain(1)_mem(2-5)-reac(6).mat'),'chain_6','-v7.3')


%% find FCSP

good=[];
for i=1:length(sessionIdx)
    chain_temp=chain_5(floor(chain_5(:,1)/100000)==sessionIdx(i),:);
    load(fullfile(homedir,'xcorr','coding_10ms','ext',sprintf('fc_decoding_f%d.mat',sessionIdx(i))));
    conn_id=cell2mat(sums0(:,2));
    %%
    good=[good;FCSP(chain_temp(:,1:2),[13,14],sums0,conn_id)&FCSP(chain_temp(:,2:3),[13,14],sums0,conn_id)...
        &FCSP(chain_temp(:,3:4),[13,14],sums0,conn_id)&FCSP(chain_temp(:,4:5),[13,14],sums0,conn_id)];
end
save(fullfile(homedir,'xcorr','chain_5_mem(1-4)-nonreac(5).mat'),'good','chain_5')


%% function
function [sig,pair]=load_sig_pair(Path,type)
if strcmp(type,'ext')
    fl=dir(fullfile(Path,'ext_conn_w_reg_*.mat'));
elseif strcmp(type,'inh')
    fl=dir(fullfile(Path,'inh_conn_w_reg_*.mat'));
end

sig=struct();
[sig.suid,sig.reg,sig.sess,sig.mem_type]=deal(cell(0));
pair=sig;
for fidx=1:length(fl)
    fstr=load(fullfile(fl(fidx,1).folder,fl(fidx,1).name)); 
    fields={'suid','reg','mem_type'};
    for fi=fields
        sig.(fi{1}){fidx}=fstr.sig_meta.(fi{1});
        pair.(fi{1}){fidx}=fstr.pair_meta.(fi{1});
    end
    sig.sess{fidx}=repmat(9+fidx,size(fstr.sig_meta.suid,1),1);
    pair.sess{fidx}=repmat(9+fidx,size(fstr.pair_meta.suid,1),1);
end        

for fi=[fields,{'sess'}]
    sig.(fi{1})=cell2mat(sig.(fi{1})');
    pair.(fi{1})=cell2mat(pair.(fi{1})'); 
end
end

function out=FCSP(in,trialIdx,sums0,conn_id)
for ii=1:size(in,1)
    sums=sums0{conn_id(:,1)==in(ii,1)&conn_id(:,2)==in(ii,2),3}(trialIdx); %[distractorNo-S1,distractorNo-S2,distractor-S1,distractor-S2]
    FCSP=cell2mat(sums);
    if sum(FCSP(1,:,12:13),'all')>0
        out(ii,:)=1;
    else
        out(ii,:)=0;
    end
end
end