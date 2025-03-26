function [series series_pred pred_err fGns fBms] = converging_prediction_log2(beta, sigma2, k, nSeries, range_w, range_delta, end_w, end_delta, min_v, doPlot, report)
% [series series_pred pred_err fGns fBms] = converging_prediction_log(beta, sigma2, k, doPlot, nSeries, range_w, range_delta, end_w, end_delta, min_v, doPlot, report)
% 
% generate fGn / fBm series that converge to a similar final value, and
% calculate the prediction for the upcoming value.
%
% converging_prediction_minsc differs from converging_prediction by
% allowing the minimum value of the series of interest to be set as an
% input (the series is linearly shifted to ensure the minimum series value
% is v_min)
%
% inputs
% ------
% beta     - beta of the fGn / fBm series
% sigma2   - variance of the fGn (default value = 1) 
% k        - # of prior samples on which the prediction for the next sample
%            is made (default = 50)
% doPlot   - set to 1 to plot results, 0 to turn off plot (default = 1)
% nSeries  - # of series to generate (default = 100)
% delta    - tolerance level for final value v_f. i.e. series are only
%            returned if |series(end) - v_f| <= delta (default = 1e-3 * sqrt(sigma2))
% v_f      - value towards which the series converge, PRIOR TO scaling for the minimum
%            (default = 1.5 * sqrt(sigma2))
% v_min    - minimum value of the series of interest (default = 0)
% progress - set to 1 to get a progress report while simulations are
%            ongoing. (default = 0)
%
%
% outputs
% -------
% series      - A matrix of size nSeries x k. Each row is a fGn/fBm of length
%               k for which the final value is within delta of v_f.
% series_pred - A vector of size nSeries x 1. series_pred(i) contains the best 
%               prediction for the value following series(i,k) (found using 
%               function fgn_pred)
% pred_err    - A vector of size nSeries x 1. pred_err(i) contains the
%               prediction error for the k+1'th sample, 
%               i.e. series_pred(i) - series(i,k+1) [series(i,k+1) not
%               included in output]


%% set default input values for unspecified inputs
if ~exist('sigma2','var') || isempty(sigma2)
    sigma2 = 1;
end

if ~exist('k','var') || isempty(k)
    k = 50;
end

if ~exist('nSeries','var') || isempty(nSeries)
    nSeries = 100;
end

if ~exist('range_w','var') || isempty(range_w)
    range_w = 4;
end

if ~exist('range_delta','var') || isempty(range_delta)
    range_delta = .05;
end

% range
if ~exist('range_w','var') || isempty(range_w)
    range_w = 4;
end

if ~exist('range_delta','var') || isempty(range_delta)
    range_delta = .01;
end

% end value
if ~exist('end_w','var') || isempty(end_w)
    end_w = log(440);
end

if ~exist('end_delta','var') || isempty(end_delta)
    end_delta = .01;
end

if ~exist('min_v','var') || isempty(min_v)
    min_v = log(220);
end

if ~exist('doPlot','var') || isempty(doPlot)
    doPlot = 1;
end

if ~exist('report','var') || isempty(report)
    report = 1;
end

if ~exist('H_delta','var') || isempty(H_delta)
    H_delta = .01;
end

%% define series' mean value

mean_v = mean([min_v min_v+log(range_w)]);

%% prepare for search

% ensure synthesized series are long enough to accommodate k
n = 10;
while 2^n + 1 < 5*k
    n = n+1;
end

% find H from beta
if beta > 1
    H = (beta-1)/2;
else
    H = (beta+1)/2;
end


%% search for series that converge to v_f, within tolerance delta

% series counter. 1 <= i_s <= nSeries
i_s = 0;

if report == 2
    fprintf(1,'\n\npercent done:\n');
end

% diagnostic variables
i_synth = 0;
report_ps = 1 : nSeries;
report_ps = downsample(report_ps, length(report_ps)/10);


% while we haven't found the # of required series...
while i_s < nSeries
    
    % synthesize an fGn/fBm
    [B x] = synthfbmcircul2(2^n+1, H, sigma2);
    i_synth = i_synth+1;
    
    % select the appropriate series for analysis
    if beta > 1
        s = B;
        isfGn = 0;
    else
        s = x;
        isfGn = 1;
    end
    
    
    % loop through the series to find a segment of length k for which the
    % k'th value is within delta of v_f
    j = k+1;
    while j+1 <= length(s)
        
        % qualification criteria:
        % - f_max = range_w * f_min --> s_max - s_min = log(range_w)
        % - f_end = end_w * f_min --> s_j - s_min = log(w)
        series_maxInd = j;
        series_minInd = j - (k-1);
        s_k = s(series_minInd : series_maxInd);
        
        mean_s = mean(s_k);
        mean_adj = mean_v; % - mean_s;
        
        rangeOK = abs( range(s_k) - log(range_w) ) <= range_delta;
        endOK   = abs( s_k(end) +mean_adj - end_w ) <= end_delta;
        minOK   = abs( min(s_k) +mean_adj - min_v ) <= end_delta;
        
        allOK = 0;
        if rangeOK && endOK && minOK

            if isfGn
                s_x = s_k;
            else
                s_x = s_k(2:end) - s_k(1:end-1);
            end
            
            H_est = wavelet_est(s_x, 1);
            
            if abs( H - H_est ) <= H_delta
                allOK = 1;
            end
        end            
        
        if allOK % rangeOK && endOK && minOK
            i_s = i_s + 1;
            
            % output progress
            if report == 2
                if any(report_ps == i_s)
                    percent_done = 100*i_s/nSeries;
                    percent_q    = 100*i_s/i_synth;

                    fprintf(1,'\n\npercent of required series found: %3.1f%%\n',percent_done);
                    fprintf(1,'percent of synthesized series that qualified so far: %3.1f%%',percent_q);                    
                end
            end
            
            
            % store the last k samples of this segment
            series(i_s,:) = s_k;
            
                                                   
            if beta > 1
                % find the best k'th order prediction for the next fGn sample
                %
                % N.B. that for fBm prediction, the index to the corresponding
                % sample of the fGn is shifted down by 1, since by
                % convention of synthfbmcircul, B[n+1] = B[n] + x[n]

                fGn_maxInd = series_maxInd - 1;
                fGn_minInd = series_minInd - 1;
                xt = x(fGn_maxInd : -1 : fGn_minInd);
                
                x_pred = fgn_pred(xt, H, sigma2);
                series_pred(i_s,1) = x_pred + series(i_s,end);
                
                fBms(i_s,:) = series(i_s,:);
                fGns(i_s,:) = x(fGn_minInd : fGn_maxInd);
            else
                % find the best k'th order prediction for the next fGn sample
                
                fGn_maxInd = series_maxInd;
                fGn_minInd = series_minInd;
                
                xt = x(fGn_maxInd : -1 : fGn_minInd);
                
                x_pred = fgn_pred(xt, H, sigma2);
                series_pred(i_s,1) = x_pred;
                
                fGns(i_s,:) = series(i_s,:);
                %%% new to _minsc
%                 fBms(i_s,:) = B(fGn_minInd+1 : fGn_maxInd+1);
                %%%
            end
            
            pred_err(i_s,1) = series_pred(i_s,1) - s(j+1);
            
            
            %%% new to _minsc
            % scale the series so its minimum value is min_v
            %
            % if the series is fBm, then the corresponding fGn x is
            % unaffected by the scaling
            %
            % if the series is fGn, then the corresponding fBm is changed
            % by the scaling of the fGn... but we don't care about this
            % case for our purposes
%             min_s = min(series(i_s,:));
%             series(i_s,:) = series(i_s,:) - min_s + min_v;
%             series_pred(i_s,1) = series_pred(i_s,1) - min_s + min_v;            
%             if beta < 1
% %                 series_pred(i_s,1) = series_pred(i_s,1) - min_s + min_v;
%                 fGns(i_s,:) = fGns(i_s,:) - min_s + min_v;
%             end
            %%%
            
% % % % %             series(i_s,:) = series(i_s,:) + mean_adj;
% % % % %             series_pred(i_s,:) = series_pred(i_s,:) + mean_adj;
% % % % %             if beta < 1
% % % % %                 fGns(i_s,:) = fGns(i_s,:) + mean_adj;
% % % % %             end
            
            
            % exit the search loop
            j = Inf;
        end

        j = j+1;
   end
end

std_mean_adj = std(series_pred)
range_mean_adj = range(series_pred)

%% report
if report >= 1
    percent_q = 100*i_s/i_synth;
    fprintf(1,'\n\npercent of synthesized series that qualified = %3.1f%%\n\n',percent_q);
end

%% plot results

if doPlot
    
    figure; hold on;
%     subplot(211); hold on;
    for h = 1:nSeries
        plot(1:k, series(h,:), 'b-');
        plot(k+1, series_pred(h), 'r.');
    end
    title(['beta = ' num2str(beta)])

%     subplot(212); hold on;
%     for h = 1:nSeries
%         plot(1:k, exp(series(h,:)), 'b-');
%         plot(k+1, exp(series_pred(h)), 'r.');
%         plot([1 k], exp(min_v)*[1 1],'g-');
%         plot([1 k], exp(min_v+log(range_w))*[1 1],'g-');
%     end    

end