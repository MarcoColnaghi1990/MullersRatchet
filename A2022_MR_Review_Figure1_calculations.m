%%%. Simulates Muller's ratchet in an asexual population or in a population
%%%. that undergoes LGT at a rate l.
%%%.
%%%. Population size = N, number of genes = g; evolve for nGen generations.
%%%.
%%%. The state of the system is described by a N x g matrix whose elements
%%%. are either 0 (no mutations) or n (number of mutations).
%%%.
%%%. Strength of selection s, fitness = (1 - s)^k
%%%. Mutation rate/genome/generation: u
%%%
%%%.

%%
clf;hold on
colors=get(gca,'colororder');
%%%. System Parameters
nGen=5000;                  % Number of generations
N=10000;                      % Population size
g=100;                       % Number of genes; 1 gene =~ 1000bp
tic
nRuns = 1;                  % # of Iterations
Lvector = [0];    % eDNA length
l = 0.0;                    % LGT rate / individual / generation
u_bp= 4*10^(-7);            % Mutation rate per bp
u = 1000*u_bp;               % Mutation rate per gene
U=u*g;                       % Genome-wide mutation rate
s = 0.01;                   % Strength of seleciton agains deleterious mutations
inLoad = 0.05;                % Initial mutation load
n_individuals_LLC = zeros(1,nGen);
n_mutations_LLC = zeros(1,nGen);
mean_mut_load = zeros(1,nGen);

%%%. Simulation
disp(['n_0 = ',num2str(1/exp(U/s) *N)])
dmdt=zeros(numel(Lvector),nRuns);

for L1 = 1:numel(Lvector)
    L=Lvector(L1);
    n_fixedM=zeros(nRuns,nGen);                 % # Fixed mutations
    totFixMut=zeros(1,nRuns);
    
    for r1 = 1:nRuns
        %%%. Initialise
        rand_mat = rand(N,g);                 % generate random matrix (for mutations)
        M = double(rand_mat < inLoad);        % N x g matrix; 0 = wild-type, 1 = mutant
        
        %%%. Evolution
        for t=2:nGen
            oldMat=M;
            mload = sum(M,1)/N;                  % compute mutation load in eDNA pool
            M = offspring(M,s);                  % next generation
            M = mutate(M,U);                    % introduce new  mutations
            
            if l>0
                M = LGT2(M,l,L,oldMat);         % LGT
            end
            
            totFixMut(r1) = sum(min(M));
            
            if t == 100
                M_100 = M;disp(t)
            end
            if t == 500
                M_500 = M;disp(t)
            end
            if t == 1000
                M_1000 = M;disp(t)
            end
            if t == 1500
                M_1500 = M;disp(t)
            end
            if t == 2000
                M_2000 = M;disp(t)
            end
            
            if mod(t,100) == 0
                clf
                disp(t)
                histogram(sum(M,2),0:25,'Normalization','pdf')
                ylim([0,0.25])
                set(gca,'FontName','Lucida Bright','Fontsize',12)
                set(gca,'YMinorTick','on','XMinorTick','on')
                ylabel('Number of individuals')
                xlabel('Mutation load')
                pause(1)
                drawnow()
            end
            
            n_fixedM(r1,t) = sum(min(M));
            if n_fixedM(r1,t)>n_fixedM(r1,t-1)
                disp(['Generation ', num2str(t), '. N fixed muts = ', num2str(n_fixedM(r1,t))])
            end
            
            n_individuals_LLC(t) = sum(sum(M,2)==min(sum(M,2)));
            n_mutations_LLC(t) = min(sum(M,2));
            mean_mut_load(t) = mean(sum(M,2));
            
        end
        
        disp(['run ',num2str(r1),' of ', num2str(nRuns), '. L = ', num2str(L), '. g = ', num2str(g), '. *** ', num2str(totFixMut(r1)), ' fixed mutations ***'])
    end
    
end
toc

clf
ax1=subplot(2,3,1);hold on
ax2=subplot(2,3,2);hold on
ax3=subplot(2,3,3);hold on
ax4=subplot(2,1,2);hold on
%ax5=subplot(4,1,4);hold on

histogram(ax1,sum(M_100,2),0:15,'Normalization','pdf')
histogram(ax2,sum(M_500,2),0:15,'Normalization','pdf')
histogram(ax3,sum(M_1000,2),0:15,'Normalization','pdf')

ax1.Position = [0.1300    0.5    0.24    0.43];
ax2.Position = [0.4100    0.5    0.24    0.43];
ax3.Position = [0.69    0.5    0.24    0.43];
ax4.Position = [0.1300    0.1100    0.7750    0.25];

for x = [ax1, ax2, ax3]
    ylim(x,[0,0.25])
    set(x,'FontName','Lucida Bright','Fontsize',12)
    set(x,'YMinorTick','on','XMinorTick','on')
end

ylabel(ax1,'Frequency');
xlabel(ax2,'Mutation load')

for x = [ax2,ax3]
    set(x,'Yticklabel',[])
end

title(ax1,'t = 100'); title(ax2,'t = 500'); title(ax3,'t = 1000');

plot(ax4,1:nGen,n_individuals_LLC,'color','k')

ylabel(ax4,'Number of individuals');ylim([0,1000])
yyaxis right
plot(ax4,1:nGen,n_mutations_LLC,'linewidth',2)

ylabel(ax4,'Mutation load');ylim([0,10])
set(ax4,'FontName','Lucida Bright','Fontsize',12)
set(ax4,'YMinorTick','on','XMinorTick','on')
xlabel(ax4,'Time (generations)')


% Create textbox
annotation('textbox',...
    [0.13 0.94 0.40 0.05],...
    'String','A',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Lucida Bright',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.13 0.36 0.40 0.05],...
    'String','B',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Lucida Bright',...
    'FitBoxToText','off');



%%%. Additional functions

function X = offspring(X,s)
N = numel(X(:,1));                            % population size
nMut = sum(X,2);                              % individual mutation load
sPower = (1-s).^nMut;
mFitness = mean(sPower);
fitVec = sPower/mFitness;                     % relative fitness
y= randsample(N,N,true,fitVec);
X = X(y,:);
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

nLGT = random('Bino',N,l);                    % # of LGT transfer events in the population
whoLGT = randsample(1:N,nLGT);                 % who undergoes LGT

if sum(nLGT)>0
    LGTloci = randi([1 g-L+1],1,sum(nLGT));   % locus where rec begins
    X(whoLGT, LGTloci:LGTloci+L-1) = oldMat(whoLGT, LGTloci:LGTloci+L-1);
    % recombination
end
end
