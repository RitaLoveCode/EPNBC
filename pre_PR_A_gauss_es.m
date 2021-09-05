%手动record:蛋白质规律
%c1:编号	c2:Y/N   c8:orth
%计算每一段的系数pp_coe/pd_coe/dp_coe/dd_coe
%-------------------------n*n------n*m-----m*n-----m*m------------------------------------------------------------
%2.（beita_1*）蛋白质得到来自蛋白质的：（根据蛋白质之间的权重分）分母：相连接的蛋白质的权重和（ssij_1行和）分子：两者之间边的权重
%3.（beita_2*）蛋白质得到来自蛋白域的：根据蛋白质(0.15)5+6（已知度+贝叶斯度）+（小贝和）/按照比例分配
%4.（beita_1*）蛋白域得到来自蛋白质的：分母：蛋白域权重和 分子：蛋白域权重
%5.（beita_2*）蛋白域得到来自蛋白域的：根据蛋白域之间的权重分
%1--------------分母：相连接的蛋白质的权重和（ssij行和）分子：两者之间边的权重----------------------
beita_1=0.75;%改动0.85-> 0.75->0.65
beita_2=0.85;
pp_coe=zeros(size(ppi));%p to p
for i=1:num_pro
    %1.
    t=find(ssij(i,:)>0);
    %2. t=find(ppi(i,:)>0);%----------------------》有改动余地（给哪些蛋白质分配）
    for j=1:length(t)
        %!!!!!!!!!!!!!!!!!!!!!!! 待改动：pp_coe(i,t(j))=ssij(i,t(j))*gauss_pp(i,t(j))/sum(gauss_pp(:,t(j)));!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %3.pp_coe(i,t(j))=ssij(i,t(j))*gauss_pp(i,t(j))/sum(gauss_pp(:,t(j)));
        %4.
        pp_coe(i,t(j))=gauss_pp(i,t(j))/sum(gauss_pp(:,t(j)));
    end
end%any(isnan([1 nan]))
pp_coe_1=beita_1*pp_coe;%beita_1
%2---------------分母：蛋白域权重和 分子：蛋白域权重-------------------------------------------------------
pd_coe=zeros(num_pro,size(dom_unique,1));%p to d
for i=1:num_pro
    t=A_pro_dom(A_pro_dom(:,1)==i,2);%蛋白质所在域
    if length(t)>0%如果属于某域：因为有的蛋白质可能不属于任何域
        for j=1:length(t)%分配原则：如果是关键蛋白，就给这个域再多分点；如果域权重和大，也多分
            pd_coe(i,t(j))=sum(wddi(t(j),:))/sum(sum(wddi(:)));
        end
    end
end%any(isnan([1 nan]))  有NAN
pd_coe_1=(1-beita_1)*pd_coe;%beita_1
%3----------------分母：蛋白域中蛋白质权重总和分子：单个蛋白质的权重和------------------------------------
dp_coe=zeros(size(dom_unique,1),num_pro);%d to p(只分配给自己域的蛋白质)
for i=1:length(dom_unique)
    t=A_pro_dom(A_pro_dom(:,2)==i,1);%d
    fenmu_3=0;
    for j=1:length(t)
        fenmu_3=fenmu_3+sum(gauss_pp(t(j),:));%分母：denominator
    end
    for j=1:length(t)
        dp_coe(i,t(j))=sum(gauss_pp(t(j),:))/fenmu_3;
    end
end%any(isnan([1 nan]))
dp_coe_1=beita_2*dp_coe;%beita_2
%4----------------蛋白域得到来自蛋白域的：根据蛋白域之间的权重分----------------------------------------
dd_coe=zeros(size(dom_unique,1),size(dom_unique,1));%d to d
for i=1:length(dom_unique)
    t=find(wddi(i,:)>0);
    for j=1:length(t)
        dd_coe(i,t(j))=wddi(i,t(j))/sum(wddi(:,t(j)));
    end
end%any(isnan([1 nan]))
dd_coe_1=(1-beita_2)*dd_coe;%beita_2
%------------------------------------------------------------------------------------------------------------------------------
%列：cat(1)  行：cat(2)
COE_1=cat(1,pp_coe,dp_coe);%列
COE_2=cat(1,pd_coe,dd_coe);%列
COE=cat(2,COE_1,COE_2);%转移矩阵大功告成

