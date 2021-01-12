function CRNSSA = CRN(C,B,L,M,N,pop_size,max_iter,dim,bound,r_a,p_c,p_m,Cmax)
%======================= INITIALIZATION===================================%
%======================= LOAD FILE LBC.txt  ====================================%;
%%%%%%%%%%%%%%%%%%%%%%%%
 B=csvread('M_B.txt');
 L=csvread('M_L.txt');
 C=csvread('M_C.txt');

 %Cmax=20;%csvread('Cmax.txt');
 [N M]=size(L); % (M): number of channel;(N): number of secondary user
C=reshape(C,N,N,M);


 pop_size=20; %number of individl
 dim=M*N; %number of dimenssion in SSA
 bound = 100; % bound of spider web
max_iter=500;  % max iteration
r_a = 1;   % rate of vibration attenuation when propagating over web
p_c = 0.7; % probability 1-p_c to change its mask
p_m = 0.1; % if mask changed, each bit of the vector has a probability of p_m to assigned=1
info = true;
maxNumRun = 30;
maxCmax = 20;
fitness_num = 1;%1: MSUM; 2:MMIN; 3:MProFair

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KQ_ssa = [];
if (info)
        %     fprintf('               SSA starts at %s\n', datestr(now));
        %     fprintf('==============================================================\n');
            
        if(fitness_num == 1)
            fprintf('                                MSUMR SSA\n');
        else if(fitness_num == 2)
            fprintf('                                MMINR SSA\n');
        else
            fprintf('                                MPFAIR SSA\n');
            end
        disp([' Num.Run=' num2str(maxNumRun)  ' Max           Mean          Std']);
    end
        fprintf('==============================================================\n');
        start_time = java.lang.System.currentTimeMillis;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Cmax=1:maxCmax
    num_KQ = [];
    for numRun = 1:maxNumRun
    %======================= CREATE x ====================================%
    x=rand(pop_size,dim)* 2 * bound - bound; %creat random matrix x(real)
    %x
    %Initialization y nhi phan theo x
    for i=1:pop_size
        for j=1:dim
            s=1/(1+exp(-x(i,j))); %S1 transfer function
            if rand < s
                y(i,j)=0;
            else
                y(i,j)=1;
            end
        end
    end

    %y(binary)=x(real)
    %======================= matxl(x,L,C,Cmax) ============================%
    y=matxl(y,L,dim,pop_size);%ktra y voi L(n,m)
    y=matxc(y,C,M,N,dim,pop_size); %ktra y voi C(n,i,m)
    y=matxcmax(y,Cmax,M,N,dim,pop_size); %ktra voi Cmax

    %======================= 12/04/2019 =================================%
    for i=1:pop_size
        for j=1:dim
            if y(i,j)==0
                y(i,j)=0;
            else
                y(i,j)=x(i,j);
            end
        end
    end

    %======================= KHOI TAO SSA =================================%

    g_best = Inf;
    g_best_hist = [];
    g_best_pos = zeros(1, dim);


    %==================== SSA TINH THEO MATRIX x ============================%
    x=y;
    target_x = x;
    target_intensity = zeros(pop_size, 1);
    mask = zeros(pop_size, dim);
    movement = zeros(pop_size, dim);
    inactive = zeros(pop_size, 1);
    start_time = java.lang.System.currentTimeMillis;
    

    iter = 0;
    while (iter < max_iter)
        iter = iter + 1;
        
        if(fitness_num == 1)
                spider_fitness = MSUMR(y,M,B,pop_size);
            else if(fitness_num == 2)
                    spider_fitness=MMINR(y,M,B,pop_size);
                else
                    spider_fitness=MPFAIR(y,N,M,B,pop_size);
                end
        end

        base_distance = mean(std(x));
        distance = squareform(pdist(x));

        if (min(spider_fitness) < g_best);
            [g_best, ind] = min(spider_fitness);
            g_best_pos = y(ind, :);   %nghiem tot nhat tim duoc

        end

        g_best_hist = [g_best_hist g_best];
%         if (info )&& ((iter == 1 || iter == 10 || (iter < 1001 && mod(iter, 100) == 0) || (iter < 10001 && mod(iter, 1000) == 0)|| (iter < 100000 && mod(iter, 1000) == 0)))
%             elapsed_time = java.lang.System.currentTimeMillis - start_time;
%             fprintf(['% 5s %.4e %.4e %.4e %.4e %02d:%02d:%02d.%03d\n'], num2str(iter), ...
%                 g_best, min(spider_fitness), base_distance, mean(mean(distance)), ...
%                 fix(elapsed_time / 3600000), mod(fix(elapsed_time / 60000), 60), ...
%                 mod(fix(elapsed_time / 1000), 60), (mod(elapsed_time, 1000)));
%         end

        intensity_source = log(1 ./ (spider_fitness + 1E-100) + 1); %fomular (1)
        intensity_attenuation = exp(-distance / (base_distance * r_a)); % formula 3
        intensity_receive = repmat(intensity_source, 1, pop_size) .* intensity_attenuation; %???

        [~, max_index] = max(intensity_receive, [], 2); %vi tri co gia tri lon nhat
        best_receive = intensity_receive(sub2ind([pop_size, pop_size], 1:pop_size, max_index'))';

        keep_target = best_receive <= target_intensity;

        keep_target_matrix = repmat(keep_target, 1, dim);
        inactive = inactive .* keep_target + keep_target; %Cs
        target_intensity = target_intensity .* keep_target + best_receive .* (1 - keep_target);

        target_x = target_x .* keep_target_matrix + x(max_index, :) .* (1 - keep_target_matrix);

        rand_x = reshape(x(sub2ind([pop_size, dim], ceil(rand(1, pop_size * dim) * pop_size), ...
            reshape(repmat(1:dim, pop_size, 1), 1, pop_size * dim))), pop_size, dim);
        new_mask = ceil(rand(pop_size, dim) + rand() * p_m - 1);
        keep_mask = rand(pop_size, 1) < p_c.^inactive;
        inactive = inactive .* keep_mask;
        keep_mask_matrix = repmat(keep_mask, 1, dim);
        mask = keep_mask_matrix .* mask + (1 - keep_mask_matrix) .* new_mask;

        follow_x = mask .* rand_x + (1 - mask) .* target_x;%formula (4)

        movement = repmat(rand(pop_size, 1), 1, dim) .* movement + ...
            (follow_x - x) .* rand(pop_size, dim);

        x = x + movement; %formula (5)

        % xu ly rang buoc cua SSA (formula (6))

        %======================= matxl(x,L,C) ================================%
        %         y=fixx(x,dim,pop_size); %convert x(real) => x(0,1)
        for i=1:pop_size
            for j=1:dim
                %s=abs((2/pi)*atan((pi/2)*x(i,j)));
                s=1/(1+exp(-x(i,j))); %S1 transfer function

                if rand < s 
                    y(i,j)=0;
                else
                    y(i,j)=1;
                end
            end
        end
        y=matxl(y,L,dim,pop_size); %ktra y voi L(n,m)
        y=matxc(y,C,M,N,dim,pop_size); %ktra y voi C(n,i,m)
        y=matxcmax(y,Cmax,M,N,dim,pop_size); %ktra voi Cmax

        %=====================================================================%

    end
    num_KQ = [num_KQ g_best];
    end % end cua for num
    % Max, Mean, Std theo Num
    max_ssa = max(1000 - num_KQ);
    mean_ssa = mean(1000 - num_KQ);
    std_ssa = std(1000 - num_KQ);
    
    
    if (info)
        if(Cmax == 5 || Cmax == 10 || Cmax == 15 || Cmax == 20)
            disp(['Cmax= ' num2str(Cmax) ':::::max_ssa= ' num2str(max_ssa) ': mean_ssa = ' num2str(mean_ssa) ':std_ssa = ' num2str(std_ssa)]);
        end
    end
    KQ_ssa = [KQ_ssa mean_ssa];
end % end CMax

%=================== MA TRAN A CAN TIM CO GTRI MIN =======================%
end_time = java.lang.System.currentTimeMillis - start_time;
disp(['Time process ' num2str(end_time)]);
%KQ = 1000- g_best;
g_best_pos = y(ind, :);   %nghiem tot nhat tim duoc
A = vec2mat(g_best_pos,M); %ma tran A can tim
CRNSSA.KQ=mean(KQ_ssa);

% CRNSSA.A=A;

plot(KQ_ssa, '-o','LineWidth',2);
%scatter(KQs,'filled');
xlabel('CMax');
if(fitness_num == 1)
        csvwrite('MSUMR_ssa.txt',KQ_ssa);
        ylabel('Max Sum Reward');
    else if(fitness_num == 2)
            csvwrite('MMINR_ssa.txt',KQ_ssa);
            ylabel('Max Min Reward');
        else
            csvwrite('MPFAIR_ssa.txt',KQ_ssa);
            ylabel('Max Proportional Fair');
        end
    end
legend ('SSA')
grid on;
hold all;
end



