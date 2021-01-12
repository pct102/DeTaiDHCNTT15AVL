function f = MPFAIR(y,N,M,B,pop_size)
% b=reshape(B',1,dim);

for i=1:pop_size
    mat = vec2mat(y(i,:),M);
    sumr=(sum(mat.*B,2)+10^-6).^(1/N);%sum reward
    PF(i)=prod(sumr); %product of sum reward
end
f=1000-reshape(PF,pop_size,1);
end
