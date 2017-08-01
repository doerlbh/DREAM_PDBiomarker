

plane1 = var(acc1_05_14_20_27_54);
s_acc1_05_14_20_27_54=zeros(size,1);

parfor n = 1:size

s_acc1_05_14_20_27_54(n) = acc1_05_14_20_27_54(n,:)*plane1.'/norm(plane1,2);

end
