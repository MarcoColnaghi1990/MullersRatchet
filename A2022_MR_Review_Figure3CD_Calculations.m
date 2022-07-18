%%%. Simulates Muller's ratchet in an asexual population or in a population
%%%. that undergoes LGT at a rate l, with discrete generations and variable
%%%. population size.
%%%.
%%%. Population size = N, number of genes = g; evolve for nGen generations.
%%%.
%%%. The state of the system is described by a N x g matrix whose elements
%%%. are either 0 (no mutations) or n (number of mutations).
%%%.
%%%. Strength of selection s, fitness = (1 - s)^k
%%%. Mutation rate/genome/generation: u
%%%.

%%
%%%. System Parameters
nGen=100000;                  % Number of generations
N_0= 1000;                      % Initial population size
K_0= N_0;
g=100;                       % Number of genes; 1 gene =~ 1000bp
tic
nRuns = 25;                  % # of Iterations
s_vector = 10.^(-5:0.1:-1);    % strength of selection
l = 1;                    % LGT rate / individual / generation
u_bp= 4*10^(-7);            % Mutation rate per bp
u = 1000*u_bp;               % Mutation rate per gene
U=u*g;                       % Genome-wide mutation rate
%s = 0.01;                   % Strength of seleciton agains deleterious mutations
inLoad = 0.0;                % Initial mutation load
f_0 = 1.05;
L = 10;

%%%. Simulation
t_ext=nGen*ones(numel(s_vector),nRuns);



for s1 = 1:numel(s_vector)
    s = s_vector(s1);
    N=zeros(nRuns,nGen);
    N(:,1)=N_0;
    n_individuals_LLC = zeros(1,nGen);
    n_mutations_LLC = zeros(1,nGen);
    mean_mut_load = zeros(1,nGen);
    
    disp(['s = ',num2str(s), ', n_0 = ',num2str(1/exp(U/s) *N_0)])
    for r1 = 1:nRuns
        %%%. Initialise
        t=1;
        rand_mat = rand(N_0,g);                 % generate random matrix (for mutations)
        M = double(rand_mat < inLoad);        % N x g matrix; 0 = wild-type, 1 = mutant
        
        %%%. Evolution
        while (t<nGen && numel(M)>0)
            oldMat=M;
            t=t+1;
            M = mutate(M,U);                    % introduce new  mutations
            if l>0
                M = LGT2(M,l,L,oldMat);         % LGT
            end
            
            M = offspring(M,s,f_0,K_0);                  % next generation
            N(r1,t)=numel(M(:,1));
        end
        t_ext(s1,r1) = t;
        
        %    disp(['t_extinction = ',num2str(t),', LLC mut load = ',num2str(n_mutations_LLC(r1,t-1))])
        disp(['t_extinction = ',num2str(t)])
    end
    
    clf;hold on
    
    %plot(s_vector,mean(t_ext_2000,2),'linewidth',2,'linestyle','none','marker','o');
    plot(10.^(-5:0.1:-1),mean(t_ext_1000,2),'linewidth',2,'linestyle','none','marker','o');
    plot(s_vector,mean(t_ext,2),'linewidth',2,'linestyle','none','marker','o');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    %legend({'K = 2000','K = 1000','K = 500'})
    drawnow
    
end

%t_ext_1000_L5_l01 = t_ext;




%%%. Additional functions

function X = offspring(X,s,f_0,K_0)
nMut = sum(X,2);                                    % individual mutation load
sPower = (1-s).^nMut;                               % individual fitness
mean_fitness = mean(sPower);                        % average fitness
old_N = numel(X(:,1));                              % old population size
N = min(K_0,sum(random('Poisson',f_0*mean_fitness,1,old_N)));    % new pop size
fitVec = sPower/mean_fitness;                           % relative fitness
y= randsample(old_N,N,true,fitVec);                     % sample individuals w/ replacement
X = X(y,:);                                         % new generations
end



%%% Mutation function. It draws a random integer k from a Poisson
%%% distribution for each individual, which is the number of new mutations
%%% acquired. For each individual, it selects k loci in its genome, which
%%% acquire a new mutation each.

function X = mutate(X,U)
N = numel(X(:,1));                            % population size
g = numel(X(1,:));                            % genome size
mutString = random('Poisson',U,N,1);

nMut = sum(mutString);
if nMut>0
    yStr = randi([1,g],1,nMut); a=1;
    pStr = zeros (1,nMut);
    mut=find(mutString>0);
    for k=1:numel(mut)
        i=mut(k);
        x=repmat(i,1,mutString(i));
        pStr(a:a+mutString(i)-1)=x;
        a=a+mutString(i);
    end
    z= sum ([g.*(pStr-1);yStr]);
    X(z) = X(z)+1;
end
%  X = min(X,1);

end


%%% LGT function. It takes as input the probability of lgt l.
%%% string length L and the old matrix.
%%%
%%% The number of individuals undergoing LGT is sampled at random from a binomial distribution
%%% with N trials and probability l; the selected individuals, sampled at
%%% random from the population, undergo LGT in the following way:
%%% A random locus x between 1 and g-L+1 is selected. Then a random row from the
%%% old matrix (corresponding to a random individual from the previous generation.
%%% Elements x to x+L-1 of the new matrix become equal to those
%%% of the previous matrix.

function X = LGT2(X,l,L,oldMat)
N = numel(X(:,1));                            % population size
g = numel(X(1,:));                            % genome size
oldN = numel(oldMat(:,1));                    % previous pop size

nLGT = random('Bino',N,l);                    % # of LGT transfer events in the population
recipients = randsample(1:N,nLGT);                 % who undergoes LGT
donors = randsample(1:oldN,nLGT);

if sum(nLGT)>0
    LGTloci = randi([1 g-L+1],1,sum(nLGT));   % locus where rec begins
    X(recipients, LGTloci:LGTloci+L-1) = oldMat(donors, LGTloci:LGTloci+L-1);
    % recombination
end
end
