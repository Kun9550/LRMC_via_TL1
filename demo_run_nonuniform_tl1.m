
    rng('default');
    set = 0;
    eg = 2;   % For uniform setting, change to eg = 1
    maxiter = 10;
    d1_values = [100, 150];
    r_values = [2, 5];
    SR_values = [0.1, 0.2];

    % Preallocate results matrices
    num_models = length(d1_values) * length(r_values) * length(SR_values);
    Erra = zeros(1, num_models);
    Tima = zeros(1, num_models);

    % Iterate over parameter combinations
    for d1 = d1_values
        d2 = d1;
        len = d1 * d2;
        for r = r_values
            for SR = SR_values
                set = set + 1;
                Err = zeros(1, maxiter);
                Tim = zeros(1, maxiter);

                % Generate random matrices once per (d1, r) combination
                ML = randn(d1, r);
                MR = randn(d2, r);
                M = ML * MR';

                % Choose the uniform/non-uniform sampling 
                if eg == 1
                    idx = randperm(len, round(SR * len));
                    options.lamfactor = 1e-3; %needed to be tuned
                    options.sigma = 5e-3; 
                elseif eg == 2
                    pvec = ones(d1, 1);
                    cnt = round(0.1 * d1);
                    pvec(1:cnt) = 2 * pvec(1:cnt);
                    pvec(cnt + (1:cnt)) = 4 * pvec(cnt + (1:cnt));
                    pvec = d1 * pvec / sum(pvec);

                    qvec = ones(d2, 1);
                    cnt = round(0.1 * d2);
                    qvec(1:cnt) = 2 * qvec(1:cnt);
                    qvec(cnt + (1:cnt)) = 4 * qvec(cnt + (1:cnt));
                    qvec = d2 * qvec / sum(qvec);

                    probmatrix = rand(d1, d2) .* (pvec * qvec');
                    [~, sortidx] = sort(probmatrix(:), 'descend');
                    nzidx = sortidx(1:round(SR * len));

                    options.lamfactor = 2e-4; %needed to be tuned
                    options.sigma = 5e-3;
                end

                % Main iteration loop
                for iter = 1:maxiter
                    % Generate observed matrix with noise
                    A = zeros(d1, d2);
                    Mobs = M(nzidx);
                    noiseLevel = 0;
                    randvec = randn(length(nzidx), 1);
                    randvec = randvec / norm(randvec);
                    A(nzidx) = Mobs + (noiseLevel * norm(Mobs)) * randvec;

                    % Call TL1 ADMM algorithm
                    t = cputime;
                    [M2, ~, ~] = TL1(A, options);
                    Tim(iter) = cputime - t;

                    % Calculate relative error
                    Err(iter) = norm(M2 - M, 'fro') / norm(M, 'fro');
                end

                % Store results
                Erra(set) = mean(Err);
                Tima(set) = mean(Tim);
            end
        end
    end

    % Display or save results as needed
    disp(Erra);
    disp(Tima);
