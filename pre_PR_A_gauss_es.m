%�ֶ�record:�����ʹ���
%c1:���	c2:Y/N   c8:orth
%����ÿһ�ε�ϵ��pp_coe/pd_coe/dp_coe/dd_coe
%-------------------------n*n------n*m-----m*n-----m*m------------------------------------------------------------
%2.��beita_1*�������ʵõ����Ե����ʵģ������ݵ�����֮���Ȩ�ط֣���ĸ�������ӵĵ����ʵ�Ȩ�غͣ�ssij_1�кͣ����ӣ�����֮��ߵ�Ȩ��
%3.��beita_2*�������ʵõ����Ե�����ģ����ݵ�����(0.15)5+6����֪��+��Ҷ˹�ȣ�+��С���ͣ�/���ձ�������
%4.��beita_1*��������õ����Ե����ʵģ���ĸ��������Ȩ�غ� ���ӣ�������Ȩ��
%5.��beita_2*��������õ����Ե�����ģ����ݵ�����֮���Ȩ�ط�
%1--------------��ĸ�������ӵĵ����ʵ�Ȩ�غͣ�ssij�кͣ����ӣ�����֮��ߵ�Ȩ��----------------------
beita_1=0.75;%�Ķ�0.85-> 0.75->0.65
beita_2=0.85;
pp_coe=zeros(size(ppi));%p to p
for i=1:num_pro
    %1.
    t=find(ssij(i,:)>0);
    %2. t=find(ppi(i,:)>0);%----------------------���иĶ���أ�����Щ�����ʷ��䣩
    for j=1:length(t)
        %!!!!!!!!!!!!!!!!!!!!!!! ���Ķ���pp_coe(i,t(j))=ssij(i,t(j))*gauss_pp(i,t(j))/sum(gauss_pp(:,t(j)));!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %3.pp_coe(i,t(j))=ssij(i,t(j))*gauss_pp(i,t(j))/sum(gauss_pp(:,t(j)));
        %4.
        pp_coe(i,t(j))=gauss_pp(i,t(j))/sum(gauss_pp(:,t(j)));
    end
end%any(isnan([1 nan]))
pp_coe_1=beita_1*pp_coe;%beita_1
%2---------------��ĸ��������Ȩ�غ� ���ӣ�������Ȩ��-------------------------------------------------------
pd_coe=zeros(num_pro,size(dom_unique,1));%p to d
for i=1:num_pro
    t=A_pro_dom(A_pro_dom(:,1)==i,2);%������������
    if length(t)>0%�������ĳ����Ϊ�еĵ����ʿ��ܲ������κ���
        for j=1:length(t)%����ԭ������ǹؼ����ף��͸�������ٶ�ֵ㣻�����Ȩ�غʹ�Ҳ���
            pd_coe(i,t(j))=sum(wddi(t(j),:))/sum(sum(wddi(:)));
        end
    end
end%any(isnan([1 nan]))  ��NAN
pd_coe_1=(1-beita_1)*pd_coe;%beita_1
%3----------------��ĸ���������е�����Ȩ���ܺͷ��ӣ����������ʵ�Ȩ�غ�------------------------------------
dp_coe=zeros(size(dom_unique,1),num_pro);%d to p(ֻ������Լ���ĵ�����)
for i=1:length(dom_unique)
    t=A_pro_dom(A_pro_dom(:,2)==i,1);%d
    fenmu_3=0;
    for j=1:length(t)
        fenmu_3=fenmu_3+sum(gauss_pp(t(j),:));%��ĸ��denominator
    end
    for j=1:length(t)
        dp_coe(i,t(j))=sum(gauss_pp(t(j),:))/fenmu_3;
    end
end%any(isnan([1 nan]))
dp_coe_1=beita_2*dp_coe;%beita_2
%4----------------������õ����Ե�����ģ����ݵ�����֮���Ȩ�ط�----------------------------------------
dd_coe=zeros(size(dom_unique,1),size(dom_unique,1));%d to d
for i=1:length(dom_unique)
    t=find(wddi(i,:)>0);
    for j=1:length(t)
        dd_coe(i,t(j))=wddi(i,t(j))/sum(wddi(:,t(j)));
    end
end%any(isnan([1 nan]))
dd_coe_1=(1-beita_2)*dd_coe;%beita_2
%------------------------------------------------------------------------------------------------------------------------------
%�У�cat(1)  �У�cat(2)
COE_1=cat(1,pp_coe,dp_coe);%��
COE_2=cat(1,pd_coe,dd_coe);%��
COE=cat(2,COE_1,COE_2);%ת�ƾ���󹦸��

