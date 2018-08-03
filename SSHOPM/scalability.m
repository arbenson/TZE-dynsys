addpath('/Users/arb/codes/TZE-dynsys/SSHOPM/tensor_toolbox-master');
basertrials = 50;
tol = 1e-6;
shift = 1;

%%
order = 3;

for dim = 5:15
    T = tenzeros([dim dim dim]);
    for i = 1:dim
        for j = (i+1):dim
            for k = (j+1):dim
                T(perms([i j k])) = (-1)^i/i + (-1)^j/j + (-1)^k/k;
            end
        end
    end

    numtrials = basetrials * 2 * dim;  % to match with dynamical systems
    evals = zeros(numtrials,1);
    successes = zeros(numtrials,1);
    begintime = tic;
    for iter = 1:numtrials
        [eval, x, flag, it] = ...
            eig_sshopm(T, 'Shift', shift, 'tol', tol);
        evals(iter) = eval;
        successes(iter) = (flag == 0);
    end
    time = toc(begintime);
    evals = evals(successes == 1);
    save(sprintf('results/SSHOPM-evals-%d-%d.mat', order, dim), ...
         'evals', 'time', 'dim', 'order');
end

%%
order = 4;

for dim = 5:15
    T = tenzeros([dim dim dim dim]);
    for i = 1:dim
        for j = (i+1):dim
            for k = (j+1):dim
                for l = (k+1):dim
                    T(perms([i j k l])) = (-1)^i/i + (-1)^j/j + (-1)^k/k + (-1)^l/l;      
                end
            end
        end
    end

    numtrials = basetrials * 2 * dim;  % to match with dynamical systems		
    evals = zeros(numtrials,1);
    successes = zeros(numtrials,1);
    begintime = tic;
    for iter = 1:numtrials
        [eval, x, flag, it] = ...
            eig_sshopm(T, 'Shift', shift, 'tol', tol);
        evals(iter) = eval;
        successes(iter) = (flag == 0);
    end
    time = toc(begintime);
    evals = evals(successes == 1);
    save(sprintf('results/SSHOPM-evals-%d-%d.mat', order, dim), ...
         'evals', 'time', 'dim', 'order');
end

%%
order = 5;

for dim = 5:15
    T = tenzeros([dim dim dim dim dim]);
    for i = 1:dim
        for j = (i+1):dim
            for k = (j+1):dim
                for l = (k+1):dim
                    for r = (l+1):dim
                        T(perms([i j k l r])) = (-1)^i/i + (-1)^j/j + (-1)^k/k + (-1)^l/l + (-1)^r/r;
                    end
                end
            end
        end
    end

    numtrials = basetrials * 2 * dim;  % to match with dynamical systems		    
    evals = zeros(numtrials,1);
    successes = zeros(numtrials,1);
    begintime = tic;
    for iter = 1:numtrials
        [eval, x, flag, it] = ...
            eig_sshopm(T, 'Shift', shift, 'tol', tol);
        evals(iter) = eval;
        successes(iter) = (flag == 0);
    end
    time = toc(begintime);
    evals = evals(successes == 1);
    save(sprintf('results/SSHOPM-evals-%d-%d.mat', order, dim), ...
         'evals', 'time', 'dim', 'order');
end
