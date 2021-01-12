function xl = matxl(xl,L,dim,pop_size)
l=reshape(L',1,dim); 
for i=1:pop_size
    for t=1:dim
        if l(t)==0
            xl(i,t)=0;
        end
    end
end
end
