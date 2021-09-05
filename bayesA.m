%ppi中无连接的，用贝叶斯预测潜在连接可能性。Pi Pj是否有连接要找到Pi Pj的共同邻居
%----------------------------------------------------------------------------
s=zeros(size(ppi));
sum_inter=sum(sum(ppi))/2;
fm=sum_inter/(num_pro*(num_pro-1)/2);
fm1=1/fm;
for i=1:num_pro
    com_i=find(ppi(i,:)==1);
    for j=1:num_pro
        tmp=1;
        if(i~=j)&&(ppi(i,j)==0)
            com_j=find(ppi(j,:)==1);
            com=intersect(com_i,com_j);%Pi Pj的共同邻居位置（即编号）
            if(length(com)>0)%可能不存在共同邻居
                for t=1:length(com)
                    Nm=0;
                    %此蛋白质所连接的蛋白质都有谁
                    link_p1=find(ppi(com(1,t),:)==1);
                    %小循环：
                    for x=1:length(link_p1)-1
                        for y=x+1:length(link_p1)
                            if ppi(link_p1(1,x),link_p1(1,y))==1
                                Nm=Nm+1;%小团体之间的已知连接
                            end
                        end
                    end%小循环#
                    Nm_=(length(link_p1)*(length(link_p1)-1) /2)- Nm;
                    tmp=tmp*fm1*((Nm+1)/(Nm_+1));
                end
                s(i,j)=fm*tmp;
            else
                s(i,j)=0;
            end
        end
    end
end%bayes：s
%-------------------------------------------------------------------------------------
sij=zeros(size(ppi));
for i=1:num_pro%归一化log x
    for j=1:num_pro
        if s(i,j)>1
            sij(i,j)=log(s(i,j));
        end
    end
end
M=max(max(sij));% M=204.7885
ssij=sij/M;
ssij(ppi==1)=1;%归一化后：log(s(i,j))/M