%ppi�������ӵģ��ñ�Ҷ˹Ԥ��Ǳ�����ӿ����ԡ�Pi Pj�Ƿ�������Ҫ�ҵ�Pi Pj�Ĺ�ͬ�ھ�
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
            com=intersect(com_i,com_j);%Pi Pj�Ĺ�ͬ�ھ�λ�ã�����ţ�
            if(length(com)>0)%���ܲ����ڹ�ͬ�ھ�
                for t=1:length(com)
                    Nm=0;
                    %�˵����������ӵĵ����ʶ���˭
                    link_p1=find(ppi(com(1,t),:)==1);
                    %Сѭ����
                    for x=1:length(link_p1)-1
                        for y=x+1:length(link_p1)
                            if ppi(link_p1(1,x),link_p1(1,y))==1
                                Nm=Nm+1;%С����֮�����֪����
                            end
                        end
                    end%Сѭ��#
                    Nm_=(length(link_p1)*(length(link_p1)-1) /2)- Nm;
                    tmp=tmp*fm1*((Nm+1)/(Nm_+1));
                end
                s(i,j)=fm*tmp;
            else
                s(i,j)=0;
            end
        end
    end
end%bayes��s
%-------------------------------------------------------------------------------------
sij=zeros(size(ppi));
for i=1:num_pro%��һ��log x
    for j=1:num_pro
        if s(i,j)>1
            sij(i,j)=log(s(i,j));
        end
    end
end
M=max(max(sij));% M=204.7885
ssij=sij/M;
ssij(ppi==1)=1;%��һ����log(s(i,j))/M