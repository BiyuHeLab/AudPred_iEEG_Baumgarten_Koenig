% for each beta,
% find the sigma^2 that yields series with range ~= log(880) - log(220) = 1.3863

% sigma2_for_H_new is distinguished from sigma2_for_H
% in that the new function creates series of length k directly from synthfbmcircul2
% (the older function generated series of length [sampling frequency *
% sequence duration] and then downsampled to k samples...)


% sigmas giving mean(range) ~= log(880) - log(220) for nTonesInSeries = 33
% beta = 0    --> sigma^2 = .11
% beta = 0.5  --> sigma^2 = .13
% beta = 1.01 --> sigma^2 = .21
% beta = 1.5  --> sigma^2 = .09
% beta = 2    --> sigma^2 = .03


% sigmas giving mean(range) ~= log(880) - log(220) for nTonesInSeries = 100
% beta = 0    --> sigma^2 = .077   
%                (with nReps = 1e6, error for mean = 0.0054735; error for median = -0.0078522)
% beta = 0.5  --> sigma^2 = .086
%                (with nReps = 1e6, error for mean = 0.0054029; error for median = -0.0082461)
% beta = 1.01 --> sigma^2 = .149
%                (with nReps = 1e6, error for mean = 0.0043193; error for median = -0.0089352)
% beta = 1.5  --> sigma^2 = .039
%                (with nReps = 1e6, error for mean = 0.034069; error for median = -0.0012529)
% beta = 2    --> sigma^2 = .00985
%                (with nReps = 1e6, error for mean = 0.078121; error for median = -0.0023737)


function sigma2_for_H_new(sigma2, beta, k, nReps)

%% beta --> H
if beta > 1
    H = (beta-1)/2;
else
    H = (beta+1)/2;
end

n = nextpow2(k);

r_d = log(880) - log(220);

r = zeros(1, nReps);
for i_rep = 1:nReps
%     i_rep
    
    [B x] = synthfbmcircul2(2^n+1, H, sigma2);

    x = x(1:k);
    B = B(1:k);
    
    if beta < 1
        s = x;
    else
        s = B;
    end
    
    r(i_rep) = range(s);
    
%     alpha = DFA_copy_edit(s,0);
%     beta_DFA(i_rep) = 2*alpha-1;
    
end


% figure; hold on;
% plot(x,'bo-');
% plot(B,'ro-');
% title(['beta = ' num2str(beta) ', DFA beta = ' num2str(beta_DFA)])

% hist(beta_DFA)

disp(' ')
disp('--------------------------')
disp(['sigma2        = ' num2str(sigma2)])
disp(['desired range = ' num2str(r_d)])
disp(' ')
disp(['mean range    = ' num2str(mean(r))])
disp(['range diff    = ' num2str(mean(r) - r_d)])
disp(' ')
disp(['median range  = ' num2str(median(r))])
disp(['range diff    = ' num2str(median(r) - r_d)])
disp('--------------------------')
disp(' ')
