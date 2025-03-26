function series = discretize(inputseries, logToneFreqs)

for i_t = 1:length(inputseries)
    [~, ind]  = min(abs(inputseries(i_t) - logToneFreqs)); %get the index that yields min value for the tone in the input series and the discretized log tone freq vector
    series(i_t) = logToneFreqs(ind);
end

end
