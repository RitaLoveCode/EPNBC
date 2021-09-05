%exp_1=zeros(num_pro,size(expression,2)-1);
%for i=1:num_pro
%for j=1:length(expression)
%       if strcmp(pros_unique{i,1},expression{j,1})
%           exp_1(i,1:36)=cell2mat(expression(j,2:37));
%        end
%    end
%end%exp_1(有exp==0的)
%---------------------------------------------------------------------------------------------------
%calculate gamad for Gaussian kernel calculation
sd=zeros(num_pro,1);
gauss_pp_e=zeros(num_pro,num_pro);
gamadd=1;
for i=1:num_pro
    sd(i,1)=norm( exp_1(i,:))^2;%norm(vector):向量的2范数
end
gamad=num_pro/sum(sd)*gamadd;%gamadd=1
%calculate Gaussian kernel for the similarity between protein: kd
for i=1:num_pro-1
    for j=i:num_pro
            gauss_pp_e(i,j)=exp(-gamad*(norm(exp_1(i,:)-exp_1(j,:)))^2);
            gauss_pp_e(j,i)=gauss_pp_e(i,j);
    end
end%gauss_pp
