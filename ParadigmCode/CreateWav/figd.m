function handle = figd(fontSize, lineWidth, markerSize)
% handle = figd(fontSize, lineWidth, markerSize)
%
% - creates a figure with default fontSize, lineWidth, and markerSize
% - if not specified, these are initialized to 15, 2, and 10 respectively


if ~exist('fontSize','var') || isempty(fontSize)
    fontSize = 15;
end

if ~exist('lineWidth','var') || isempty(lineWidth)
    lineWidth = 2;
end

if ~exist('markerSize','var') || isempty(markerSize)
    markerSize = 10;
end

handle = figure;
set(handle, 'DefaultAxesFontSize', fontSize);
set(handle, 'DefaultLineLineWidth', lineWidth);
set(handle, 'DefaultLineMarkerSize', markerSize);