function [a_output, b_output] = bisection_function_vertical(a_input, b_input)

global Tolerance_vertical h adjacency_matrix partition manifolds n

a = a_input;
b = b_input;

while abs(b-a) > Tolerance_vertical
% bisection method to find the upper bound of the synchronization area, take value in the center of
% a-b and then determine whether it synchronizes. Determine new local bounds. Iterate until the
% distance between a and b is less than the specified Tolerance.

    h = (a+b)/2;

    synch = synchronization_function(adjacency_matrix,manifolds,partition,n);

    % determine new local bounds
    if synch
    % if it synchronizes, the centerpoint is the new local lower bound
        a = h;
    else
    % else the centerpoint is the new upper bound
        b = h;
    end
end

a_output = a;
b_output = b;

end