
%function MLBC=LBC(N,M,K,dmax,dmin,DPR)
%======================= INITIALIZATION===================================%
clc;
clear;
close all;
M=20; % number of channel
N=20; % number of secondary user
K=20; % number of primary user
dmax=4; % the maximum transmission power of secondary users
dmin=1; % the minimum transmission power of secondary users
DPR =2; % interference range (protection area)
MaxX = 10;
MaxY = 10;

sumAll = 0;

while sumAll < 1
    mang = 1:M;
    for kpu=1:K % PU k
        xk(kpu)= randi([0 MaxX]);
        index = randi([1 numel(mang)]);
        mk(kpu) = mang(index);
        mang(index) = [];
        canDuoi = uint8((MaxX/M)*(mk(kpu) - 1));
        canTren = uint8((MaxX/M)*(mk(kpu)));
        %canDuoi = (mk(kpu) - 1);
        %canTren = (mk(kpu));
        yk(kpu)= randi([canDuoi canTren]);

        % Khoi tao mang toa do nguoi dung su n tren kenh m N(n,2,m)
        for nsu = 1:N
            % kenh m la max
            if(canTren == MaxY)
                NC(nsu, 1,mk(kpu)) = randi([0 MaxX]) + rand; 
                NC(nsu, 2,mk(kpu)) = randi([0 canDuoi])+ rand;
            % kenh m la min
            else if(canDuoi == 0)
                    NC(nsu, 1,mk(kpu)) = randi([0 MaxX])+ rand; 
                    if(canTren > MaxY)
                        NC(nsu, 2,mk(kpu)) = MaxY;
                    else
                        NC(nsu, 2,mk(kpu)) = randi([canTren MaxY])+ rand;
                    end
                else
            % kenh m nam giua
                    NC(nsu, 1,mk(kpu)) = randi([0 MaxX])+ rand; 
                    if(randi([0 1]) == 1)
                        NC(nsu, 2,mk(kpu)) = randi([0 canDuoi])+ rand;
                    else
                        if(canTren > MaxY)
                            NC(nsu, 2,mk(kpu)) = MaxY;
                        else
                            NC(nsu, 2,mk(kpu)) = randi([canTren MaxY])+ rand;
                        end
                    end
                end
            end        
        end
    end

    % Ma tran khoang cach giua K primary users
    for i=1:K-1 % PU k
        for j=i+1:K
            D(i,j)=sqrt((xk(i)-xk(j))^2+(yk(i)-yk(j))^2);
            D(j,i)=D(i,j);
        end
    end



    % Ma tran khoang cach giua N secondary users
    for m = 1:M
        for i=1:N-1 % PU k
            for j=i+1:N
                DU(i,j,m)=sqrt((NC(i,1,m)-NC(j,1,m))^2+(NC(i,2,m)-NC(j,2,m))^2);
                DU(j,i,m)=DU(i,j,m);
            end
        end
    end

    for n=1:N %N dong
        for m=1:M
            mangMax = max(D);
            benhat=max(mangMax);
            % flag xet chung kenh: neu flag = 0 nghia la 
            % n voi xk nam chung 1 kenh
            flag=0;
            for k=1:K %K cot

                DIST(n,k)=0;
                % Xet xem toa do nguoi dung k va n co chung kenh ko
                bandwith = abs(NC(n, 2, m)-yk(k));
                if((bandwith == 10 && yk(k)\10 == 0) || ((NC(n, 2, m)/10) == yk(k)/10))
                    flag = 0;
                else
                    flag=1;
                    DIST(n, k)= sqrt((NC(n, 1, m)-xk(k))^2+(NC(n, 2, m)-yk(k))^2) - DPR;

                    if (benhat > DIST(n,k))
                        benhat = DIST(n,k);
                    end
                end  
            end
            if flag
                % DSE chua bang thong tren thong luong co the dat duoc
                % la duong di ngan tu nguoi dung thu n tren den kenh m
                DSE(n,m)=min(dmax,benhat);
            else
                DSE(n,m)=min(dmax,0);
            end
            %%%%%%%%%%%%%%%%%%% TINH DSEnm %%%%%%%%%%%%%%%%%%%%%%%%%

            % Gan gia tri kenh kha dung cua n tren kenh m trong ma tran L
            if DSE(n,m) > dmin
                B(n,m)=DSE(n,m).^2;
                L(n,m)= 1;
            else
                B(n,m)=0;
                L(n,m)=0;
            end
        end

    end

    try
        for n=1:N-1
            for i=n+1:N
                for m=1:M

                    DISTni(n,i,m)= DU(n,i,m);
                    % 2 kenh n co the join vo kenh chinh m
                    if  DSE(n,m)+DSE(i,m) >= DISTni(n,i,m)
                        C(n,i,m)= 1;
                        C(i,n,m)=C(n,i,m);
                    else
                        C(n,i,m)= 0;
                        C(i,n,m)=C(n,i,m);

                    end
                end
            end
        end
    catch e
        throw(e)
    end

    xn=[];
    yn=[];
    % Tao mang toa do nguoi dung n; toa do nguoi dung su n tren kenh m N(n,2,m)
    for m = 1:M
        for n=1:N
            xn = [xn NC(n, 1,m)];
            yn = [yn NC(n, 2,m)];
        end
    end

    sumAll = sum(sum(B, 'omitnan'));
    if sumAll >= 55 && sumAll <= 65
        B = B*2;
    end
    if sumAll >= 110 && sumAll <= 130
        sumRow =  sum(B');%+ sum(B');
        for i = 1:M
            if sumRow(i) == 0
                sumAll = -1;
                break;
            end
        end
    else
        sumAll
        sumAll = -1;
    end
    
end

sumAll
csvwrite('M_L.txt',L);
csvwrite('M_B.txt',B);
csvwrite('M_C.txt',C);
csvwrite('xk.txt',xk);
csvwrite('yk.txt',yk);
csvwrite('xn.txt',xn);
csvwrite('yn.txt',yn);
