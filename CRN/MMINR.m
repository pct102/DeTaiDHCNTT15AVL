function f = MMINR(y,M,B,pop_size)
% b=reshape(B',1,dim);

for i=1:pop_size
    mat = vec2mat(y(i,:),M);
    sumr=sum(mat.*B,2); %sum reward
    sr(i)=min(sumr); % min sum reward
end
f=1000-reshape(sr,pop_size,1);
end
