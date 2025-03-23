function [series_pred, w] = find_predicted_tone(series, beta)

if beta < 1
    is_fGn = 1;
    H = (beta+1)/2;

    x = log(series);    
    x = x - log(440);
else
    is_fGn = 0;
    H = (beta-1)/2;

    x = diff(log(series));
end

sigma2=1;

xrev = x(end:-1:1);
[x_pred, w, sigma2_e] = fgn_pred(xrev, H, sigma2);

if is_fGn
    series_pred = exp(x_pred + log(440));
else
    series_pred = exp(x_pred + log(series(end)));
end