%�в���beita
%A_pro_dom(1857)pro-dom
%A_dd_inter(4796)dom-dom
%ddi(744x744)
%wdd(i,a(j))=beita*cnt_link+(1-beita)*cnt_wei;%���������ӣ�������Ȩ�غͣ�
%--------------------------------------------------------------------------
%P_inD=zeros(length(dom_unique),1);
%for ii=1:length(dom_unique)
%    cnt=0;
%for jj=1:length(A_pro_dom)
%       if ii==A_pro_dom(jj,2)
%           cnt=cnt+1;
%           P_inD(ii,cnt)=A_pro_dom(jj,1);
%       end
%    end
%end%#����P_inD
%--------------------------------------------------------------------------
beita_0=0.05;
%beita_00=0.1;%�Ķ�5/10
wdd=zeros(length(dom_unique),length(dom_unique));
for i=1:length(dom_unique)
    obj1=P_inD(i,:);obj1(obj1==0)=[];%P_inD������
    dd=find(ddi(i,:)==1);%�������ӵ���֮���Ȩ��
    for j=1:length(dd)%����1���������Ȩ��
        cnt_link=0;cnt_wei=0;%����������������Ȩ�غ�
        obj2=P_inD(dd(j),:);obj2(obj2==0)=[];
        for ii=1:length(obj1)%obj=P_inDȥ��0Ԫ��֮��Ϊ���а����ĵ�����
            for jj=1:length(obj2)
                %0.if ssij(obj1(ii),obj2(jj))~=0
                    %0.
                    if ppi(obj1(ii),obj2(jj))~=0
                    %1.cnt_link=cnt_link+1;
                    %1.cnt_wei=cnt_wei+beita_00*ssij(obj1(ii),obj2(jj));
                    %!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ķ�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %2.cnt_wei=cnt_wei+ssij(obj1(ii),obj2(jj))*gauss_pp(obj1(ii),obj2(jj));
                    %!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ķ�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %3.
                    cnt_wei=cnt_wei+gauss_pp(obj1(ii),obj2(jj));%3
                    %!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ķ�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                end
                %!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ķ�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if ppi(obj1(ii),obj2(jj))==1
                    cnt_link=cnt_link+1;
                end
                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!�Ķ�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end
        end
        wdd(i,dd(j))=beita_0*cnt_link+(1-beita_0)*cnt_wei;%ddȨ�غ��Ĺ�ʽ
    end
end%#����wdd
wdd_max=max(max(wdd));%194.3048
wddi=wdd/wdd_max;%��һ��:x-min/max-min;min==0
