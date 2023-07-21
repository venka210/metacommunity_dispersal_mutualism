function [value,isterminal,direction] = nonNegativeEvent(t,y)
% Set the threshold value for non-negativity
threshold = 0;

% Determine the value of y - threshold
value = y - threshold;

% Set the terminal condition
isterminal = 1;

% Set the direction of the event
direction = -1;
end