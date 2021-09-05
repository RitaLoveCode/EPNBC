%手动对a排序之后
is_=length(find(a(1:    round(num_pro*0.01)     ,       13)==1));
p_=is_/round(num_pro*0.01);
disp(['1%个数：',num2str(is_),'     成功率：',num2str(p_)]);

is_=length(find(a(1: round(num_pro*0.05) ,13)==1));
p_=is_/round(num_pro*0.05) ;
disp(['5%个数：',num2str(is_),'     成功率：',num2str(p_)]);

is_=length(find(a(1: round(num_pro*0.1) ,13)==1));
p_=is_/round(num_pro*0.1) ;
disp(['10%个数：',num2str(is_),'     成功率：',num2str(p_)]);

is_=length(find(a(1: round(num_pro*0.15) ,13)==1));
p_=is_/round(num_pro*0.15) ;
disp(['15%个数：',num2str(is_),'     成功率：',num2str(p_)]);

is_=length(find(a(1: round(num_pro*0.2) ,13)==1));
p_=is_/round(num_pro*0.2);
disp(['20%个数：',num2str(is_),'     成功率：',num2str(p_)]);

is_=length(find(a(1: round(num_pro*0.25) ,13)==1));
p_=is_/round(num_pro*0.25);
disp(['25%个数：',num2str(is_),'     成功率：',num2str(p_)]);