function [value, isterminal, direction] = eventFunction(t, y, threshold)
    % Check if any variable is below the threshold
    belowThreshold = any(y < threshold);

    % Set the event value to the logical condition
    value = belowThreshold;

    % Terminate integration when the event is triggered
    isterminal = 1;

    % Direction doesn't matter in this case, so set it to 0
    direction = 0;
end