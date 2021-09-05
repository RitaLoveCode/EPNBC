%calculate gamad for Gaussian kernel calculation
sd=zeros(num_pro,1);
gauss_pp_s=zeros(num_pro,num_pro);
gamadd=1;
for i=1:num_pro
    sd(i,1)=norm( ssij(i,:))^2;%norm(vector):向量的2范数
end
gamad=num_pro/sum(sd)*gamadd;%gamadd=1
%calculate Gaussian kernel for the similarity between protein: kd
for i=1:num_pro-1
    for j=i:num_pro
            gauss_pp_s(i,j)=exp(-gamad*(norm(ssij(i,:)-ssij(j,:)))^2);
            gauss_pp_s(j,i)=gauss_pp_s(i,j);
    end
end%gauss_pp
