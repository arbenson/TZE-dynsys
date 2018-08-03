addpath('/Users/arb/codes/TZE-dynsys/SDP/AReigSTensors');
addpath('/Users/arb/codes/TZE-dynsys/SDP/sedumi-master/');
addpath('/Users/arb/codes/TZE-dynsys/SDP/gloptipoly3');
%%
order = 3;  % order of tensor
p = 2;      % z-eigenvector

for dim = 12:15
    mpol('x', dim);
    f_A = x(1)^0 - 1;  % start with 0 polynomial
    for i = 1:dim
        for j = 1:dim
            for k = 1:dim
                v = (-1)^i/i + (-1)^j/j + (-1)^k/k;
                f_A = f_A + v * x(i) * x(j) * x(k);
            end
        end
    end

    begintime = tic;
    f_c = [];
    profile on
    [evals, eigvec, info] = AReigSTensors(f_A, f_c, mpol(x), order, p, []);
    time = toc(begintime);
    save(sprintf('results/SDP-evals-%d-%d.mat', order, dim), ...
         'evals', 'time', 'dim', 'order');
    mset clear;
end

%%
order = 4;  % order of tensor
p = 2;      % z-eigenvector

for dim = 5:15
    mpol('x', dim);
    f_A = x(1)^0 - 1;  % start with 0 polynomial
    for i = 1:dim
        for j = 1:dim
            for k = 1:dim
                for l = 1:dim
                    v = (-1)^i/i + (-1)^j/j + (-1)^k/k + (-1)^l/l;
                    f_A = f_A + v * x(i) * x(j) * x(k) * x(l);
                end
            end
        end
    end

    begintime = tic;
    f_c = [];
    profile on
    [evals, eigvec, info] = AReigSTensors(f_A, f_c, mpol(x), order, p, []);
    time = toc(begintime);
    save(sprintf('results/SDP-evals-%d-%d.mat', order, dim), ...
         'evals', 'time', 'dim', 'order');
    mset clear;
end

%%
order = 5;  % order of tensor
p = 2;      % z-eigenvector

for dim = 5:15
    mpol('x', dim);
    f_A = x(1)^0 - 1;  % start with 0 polynomial
    for i = 1:dim
        for j = 1:dim
            for k = 1:dim
                for l = 1:dim
                    for r = 1:dim
                        v = (-1)^i/i + (-1)^j/j + (-1)^k/k + (-1)^l/l + (-1)^r/r;
                        f_A = f_A + v * x(i) * x(j) * x(k) * x(l) * x(r);
                    end
                end
            end
        end
    end

    begintime = tic;
    f_c = [];
    profile on
    [evals, eigvec, info] = AReigSTensors(f_A, f_c, mpol(x), order, p, []);
    time = toc(begintime);
    save(sprintf('results/SDP-evals-%d-%d.mat', order, dim), ...
         'evals', 'time', 'dim', 'order');
    mset clear;
end
