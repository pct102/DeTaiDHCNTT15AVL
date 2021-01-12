function x=matxcmax(x,Cmax,M,N,dim,pop_size)
for y=1:pop_size
    mat=vec2mat(x(y,:),M);
    for n=1:N
        z=sum(mat,2);
        if z(n) > Cmax
            r=zeros(1,z(n)-Cmax);%random
            for i=1:z(n)-Cmax
                flag=1;
                while flag
                    a=randi(M,1,1);
                    flag=0;
                    for j=1:z(n)-Cmax
                        if r(j)== a
                            flag=1;
                        end
                    end
                    if (mat(n,a)== 0)
                        flag=1;
                    end
                end
                r(i)=a;
            end           
            for t=1:numel(r)              
                mat(n,r(t))=0;
            end
        end    
    end %
    x(y,:)=reshape(mat',1,dim);
end
end