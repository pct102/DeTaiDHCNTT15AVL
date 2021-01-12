%{

%}
function CRNSSA = CRN(C,B,L,M,N,pop_size,max_iter,dim,bound,r_a,p_c,p_m,Cmax)
%======================= INITIALIZATION===================================%
clc;
clear;
close all;
%======================= LOAD FILE LBC.txt  ====================================%;
%%%%%%%%%%%%%%%%%%%%%%%%
 B=csvread('M_B.txt');
 L=csvread('M_L.txt');
 C=csvread('M_C.txt');

 %Cmax=20;%csvread('Cmax.txt');
 [N M]=size(L); % (M): number of channel;(N): number of secondary user
%K=numel(xk); % number of primary user
C=reshape(C,N,N,M);

fitness_num = 1;%1: MSUM; 2:MMIN; 3:MProFair
 pop_size=20; %number of individl
 dim=M*N; %number of dimenssion in ACO
maxIter = 500;  % max iteration
info = true;
maxNumRun = 30;
maxCmax = 20;

%% Problem Definition

% Create the graph 
[ graph ]  = createGraph();

start_time = java.lang.System.currentTimeMillis;
if (info)
    %     fprintf('               SSA starts at %s\n', datestr(now));
    %     fprintf('==============================================================\n');
    %fprintf(' iter g_best     pop_min    base_dist  mean_dist  time_elapsed\n');
    fprintf('                      SSA                                              ACO\n');
    disp([' Num.Run=' num2str(maxNumRun)  ' Max         Mean        Std                Max             Mean           Std']);
    fprintf('==============================================================\n');
    start_time = java.lang.System.currentTimeMillis;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACO Parameters

antNo = pop_size;
tau0 = 100 * 1 / (  graph.n * mean( graph.edges(:)  )  );  % Initial phromone concentration

eta = 1./ graph.edges;  % desirability of each edge 

rho = 0.05; % Evaporation rate 
alpha = 1;  % Phromone exponential parameters 
beta = 1;  % Desirability exponetial paramter
KQ_aco = [];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Cmax=1:maxCmax
    num_KQ_aco = [];
    for numRun = 1:maxNumRun
        %disp(['Cmax =' num2str(Cmax) '; numRun =' num2str(numRun)]);


        %% Initialization
        tau = tau0 * ones( graph.n , graph.n); % Phromone matirx 
        x = tau;
        %======================= CREATE x ====================================%
        %Initialization y nhi phan theo x cho ACO & SSA
        for i=1:pop_size
            for j=1:dim
                % acO
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

        g_best = Inf;
        g_best_hist = [];
        g_best_pos = zeros(1, dim);


    %==================== ACO TINH THEO MATRIX x ============================%
        x=y;
        bestFitness = inf;
        bestTour = [];

        iter = 0;
    % Loop CMax 1: 20

        %% ACO Main Loop
        while (iter < maxIter)
            iter = iter + 1;

            if(fitness_num == 1)
                fitness = MSUMR(y,M,B,pop_size);
            else if(itness_num == 2)
                    fitness=MMINR(y_aco,M,B,pop_size);
                else
                    fitness=MPFAIR(y_aco,N,M,B,pop_size);
                end
            end


            if (min(fitness) < g_best);
                [g_best, ind] = min(fitness);
                g_best_pos = y(ind, :);   %nghiem tot nhat tim duoc

            end

            g_best_hist = [g_best_hist g_best];
            
            % Create Ants 
        colony = [];
        colony = createColony( graph, colony , antNo, tau, eta, alpha,  beta);


        % Calculate the fitness values of all ants 
        for i = 1 : antNo 
            colony.ant(i).fitness = fitnessFunction(colony.ant(i).tour , graph );
        end

        % Find the best ant
        allAntsFitness = [ colony.ant(:).fitness ];
        [ minVal , minIndex ] = min( allAntsFitness );
        if minVal < bestFitness 
            bestFitness = colony.ant(minIndex).fitness;
            bestTour = colony.ant(minIndex).tour;
        end

        % Update phromone matrix 
        tau = updatePhromone( tau , colony );  

        % Evaporation 
        tau  = ( 1 - rho ) .* tau;
        x = tau;

        % Display the results 

        outmsg = [ 'Iteration #' , num2str(iter) , ' Shortest length = ' , num2str(bestFitness)  ];
        disp(outmsg)

            %======================= matxl(x,L,C) ================================%
            %         y=fixx(x,dim,pop_size); %convert x(real) => x(0,1)
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

            y=matxl(y,L,dim,pop_size); %ktra y voi L(n,m)
            y=matxc(y,C,M,N,dim,pop_size); %ktra y voi C(n,i,m)
            y=matxcmax(y,Cmax,M,N,dim,pop_size); %ktra voi Cmax

            %=====================================================================%

        end
        num_KQ_aco = [num_KQ_aco g_best];
    %% Results
    end % end cua for num
    % Max, Mean, Std theo Num
%     max_ssa = max(1000 - num_KQ);
%     mean_ssa = mean(1000 - num_KQ);
%     std_ssa = std(1000 - num_KQ);
    
    max_aco = max(1000 - num_KQ_aco);
    mean_aco = mean(1000 - num_KQ_aco);
    std_aco = std(1000 - num_KQ_aco);
    
    if (info)
        if(Cmax == 5 || Cmax == 10 || Cmax == 15 || Cmax == 20)
            disp(['Cmax= ' num2str(Cmax)  ':::::max_aco= ' num2str(max_aco) ': mean_aco = ' num2str(mean_aco) ':std_aco = ' num2str(std_aco)]);
        end
    end
    KQ_aco = [KQ_aco mean_aco];
end % end cua for Cmax

%=================== MA TRAN A CAN TIM CO GTRI MIN =======================%
end_time = java.lang.System.currentTimeMillis - start_time;
disp(['Time process ' num2str(end_time)]);
% KQ = 1000- g_best;
g_best_pos_aco = y(ind, :);   %nghiem tot nhat tim duoc
A = vec2mat(g_best_pos_aco,M); %ma tran A can tim
% CRNSSA.KQ=KQ;
   

% CRNSSA.A=A;
% Plot SSA
plot(KQ_aco, '-x','LineWidth',2);
hold on
 xlabel('CMax');
  if(fitness_num == 1)
        csvwrite('MSUMR_aco.txt',KQ_aco);
        ylabel('Max Sum Reward');
    else if(fitness_num == 2)
            csvwrite('MMINR_aco.txt',KQ_aco);
            ylabel('Max Min Reward');
        else
            csvwrite('MPFAIR_aco.txt',KQ_aco);
            ylabel('Max Proportional Fair');
        end
    end
%legend ('ACO','SSA', 'Location', 'southeast');
 grid on;
hold all;
hold off
end



