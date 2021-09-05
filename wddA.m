%有参数beita
%A_pro_dom(1857)pro-dom
%A_dd_inter(4796)dom-dom
%ddi(744x744)
%wdd(i,a(j))=beita*cnt_link+(1-beita)*cnt_wei;%蛋白质连接（数量和权重和）
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
%end%#产生P_inD
%--------------------------------------------------------------------------
beita_0=0.05;
%beita_00=0.1;%改动5/10
wdd=zeros(length(dom_unique),length(dom_unique));
for i=1:length(dom_unique)
    obj1=P_inD(i,:);obj1(obj1==0)=[];%P_inD行向量
    dd=find(ddi(i,:)==1);%求有连接的域之间的权重
    for j=1:length(dd)%求域1与其他域的权重
        cnt_link=0;cnt_wei=0;%连接数量、蛋白质权重和
        obj2=P_inD(dd(j),:);obj2(obj2==0)=[];
        for ii=1:length(obj1)%obj=P_inD去除0元素之后，为域中包含的蛋白质
            for jj=1:length(obj2)
                %0.if ssij(obj1(ii),obj2(jj))~=0
                    %0.
                    if ppi(obj1(ii),obj2(jj))~=0
                    %1.cnt_link=cnt_link+1;
                    %1.cnt_wei=cnt_wei+beita_00*ssij(obj1(ii),obj2(jj));
                    %!!!!!!!!!!!!!!!!!!!!!!!!!!!改动!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %2.cnt_wei=cnt_wei+ssij(obj1(ii),obj2(jj))*gauss_pp(obj1(ii),obj2(jj));
                    %!!!!!!!!!!!!!!!!!!!!!!!!!!!改动!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    %3.
                    cnt_wei=cnt_wei+gauss_pp(obj1(ii),obj2(jj));%3
                    %!!!!!!!!!!!!!!!!!!!!!!!!!!!改动!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                end
                %!!!!!!!!!!!!!!!!!!!!!!!!!!!改动!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if ppi(obj1(ii),obj2(jj))==1
                    cnt_link=cnt_link+1;
                end
                %!!!!!!!!!!!!!!!!!!!!!!!!!!!!改动!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end
        end
        wdd(i,dd(j))=beita_0*cnt_link+(1-beita_0)*cnt_wei;%dd权重核心公式
    end
end%#产生wdd
wdd_max=max(max(wdd));%194.3048
wddi=wdd/wdd_max;%归一化:x-min/max-min;min==0
