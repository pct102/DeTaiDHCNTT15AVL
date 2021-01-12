function xc=matxc(xc,C,M,N,dim,pop_size)
for i=1:pop_size
    mat=vec2mat(xc(i,:),M);
    for n=1:N
        for k=n
            for m=1:M
                if (C(n,k,m)== 1 && mat(n,m) == 1 && mat(k,m) == 1)
                    v=rand();
                    if v<0.5
                        mat(n,m)= 0;
                    else
                        mat(k,m)=0;
                    end
                end
            end
        end
    end %end for
    xc(i,:)=reshape(mat',1,dim);
end

end
