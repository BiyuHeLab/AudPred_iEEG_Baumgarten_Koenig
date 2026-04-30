%% --- Nonparametric 2-way RM ANOVA using permutation testing ---
% This function takes the same inputs as your rm_anova2 call:
% DV = GroupVec_FTPLrating{i_tone};          (vector of all ratings)
% SUB = GroupVec_sub{i_tone};                (subject ID per trial)
% A   = GroupVec_Predp34{i_tone};            (expected FTP levels 1–3)
% B   = GroupVec_p34{i_tone};                (presented FTP levels 1–4)

function stats = perm_rm_anova2(DV, SUB, A, B, n_perm)

    if nargin < 5
        n_perm = 5000;
    end

    subjects = unique(SUB);
    nSub = length(subjects);

    levelsA = unique(A);
    levelsB = unique(B);

    nA = length(levelsA);
    nB = length(levelsB);

    % --- STEP 1: convert data into subject × A × B matrix ---
    Y = nan(nSub, nA, nB);

    for i = 1:nSub
        for ia = 1:nA
            for ib = 1:nB
                filt = SUB == subjects(i) & A == levelsA(ia) & B == levelsB(ib);
                Y(i, ia, ib) = mean(DV(filt));  % collapse trials
            end
        end
    end

    % --- STEP 2: compute true F-values (subject-centered) ---
    [F_A_true, F_B_true, F_inter_true] = compute_rm_F(Y);

    % --- STEP 3: permutation distribution ---
    F_A_perm = zeros(n_perm,1);
    F_B_perm = zeros(n_perm,1);
    F_inter_perm = zeros(n_perm,1);

    for p = 1:n_perm
        % permute within subject: shuffle condition labels A,B jointly
        Yp = Y;
        for i = 1:nSub
            idx = randperm(nA*nB);
            Yp(i,:,:) = reshape(Y(i,idx), [nA nB]);
        end

        [Fa, Fb, Fi] = compute_rm_F(Yp);
        F_A_perm(p) = Fa;
        F_B_perm(p) = Fb;
        F_inter_perm(p) = Fi;
    end

    % --- STEP 4: p-values ---
    pA = mean(F_A_perm >= F_A_true);
    pB = mean(F_B_perm >= F_B_true);
    pI = mean(F_inter_perm >= F_inter_true);

    % --- return struct ---
    stats.F_A = F_A_true;
    stats.F_B = F_B_true;
    stats.F_interaction = F_inter_true;

    stats.p_A = pA;
    stats.p_B = pB;
    stats.p_interaction = pI;

    stats.levelsA = levelsA;
    stats.levelsB = levelsB;
    stats.n_perm = n_perm;
end


%% Helper: compute repeated-measures F statistics
function [F_A, F_B, F_inter] = compute_rm_F(Y)

    % Y: subjects × A × B

    [nSub, nA, nB] = size(Y);

    grandMean = mean(Y(:));

    % main effects means
    meanA = squeeze(mean(mean(Y,3),1)); % vector nA
    meanB = squeeze(mean(mean(Y,2),1)); % vector nB

    % interaction means
    meanAB = squeeze(mean(Y,1));        % nA × nB

    % subject means
    meanS = squeeze(mean(mean(Y,3),2)); % nSub × 1

    % --- Sums of squares (RM design) ---
    SSA = nSub * nB * sum((meanA - grandMean).^2);
    SSB = nSub * nA * sum((meanB - grandMean).^2);

    % interaction
    ss_inter = 0;
    for ia = 1:nA
        for ib = 1:nB
            ss_inter = ss_inter + nSub * (meanAB(ia,ib) - meanA(ia) - meanB(ib) + grandMean).^2;
        end
    end
    SSAB = ss_inter;

    % error terms
    SSE_A = 0;
    SSE_B = 0;
    SSE_AB = 0;

    for s = 1:nSub
        for ia = 1:nA
            SSE_A = SSE_A + nB*( mean(squeeze(mean(Y(s,ia,:),3))) - meanA(ia) - meanS(s) + grandMean )^2;
        end
        for ib = 1:nB
            SSE_B = SSE_B + nA*( mean(squeeze(mean(Y(s,:,ib),2))) - meanB(ib) - meanS(s) + grandMean )^2;
        end
        for ia = 1:nA
            for ib = 1:nB
                err = Y(s,ia,ib) - meanA(ia) - meanB(ib) - meanS(s) + 2*grandMean;
                SSE_AB = SSE_AB + err^2;
            end
        end
    end

    % --- F-values ---
    F_A = SSA / SSE_A;
    F_B = SSB / SSE_B;
    F_inter = SSAB / SSE_AB;
end