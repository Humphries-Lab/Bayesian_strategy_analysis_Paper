function [next_x_current,next_x_new] = abrupt_switch

% abrupt switch: immediately executes the new strategy and no other from
% trial t
next_x_current = 0;
next_x_new = 1;