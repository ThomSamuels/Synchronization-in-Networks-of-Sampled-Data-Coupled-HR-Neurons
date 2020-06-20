function [a_output, b_output] = bisection_function_horizontal(a_input, b_input)

global adjacency_matrix manifolds partition
global Tolerance_horizontal sigma n

a = a_input;
b = b_input;

while abs(b-a) > Tolerance_horizontal
% bisection method to find the left bound of the synchronization area, take value in the center of
% a-b and then determine whether it synchronizes. Determine new local bounds. Iterate until the
% distance between a and b is less than the specified Tolerance.

    sigma = (a+b)/2;

    synch = synchronization_function(adjacency_matrix,manifolds,partition,n);

    % determine new local bounds
    if synch
    % if it synchronizes, the centerpoint is the new local right bound
        b = sigma;
    else
    % else the centerpoint is the new left bound
        a = sigma;
    end
end

a_output = a;
b_output = b;

end