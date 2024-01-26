function [value,isterminal,direction] = nonNegativeEvent(t,y)
% Set the threshold value for non-negativity
threshold = 10^-6;
%y(3:3:end,:) 
z = 1-any(abs(y(3:3:end,:))<threshold); % this looks at all the 'm' values in the inputted y array in the ode solver. here they start from the 3rd column and are every third column until the end
% Determine the value of y - threshold
value = double(z);%y(3,:)-threshold %this converts the logical value of z into a double that can be evaluated

% Set the terminal condition
isterminal = 1;

% Set the direction of the event
direction = 0;
end