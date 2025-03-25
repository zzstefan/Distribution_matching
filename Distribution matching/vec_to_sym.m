function harmonized_matrix = vec_to_sym(vec, n)
    % vec is the matrix where each row is a vector to be transformed back to a symmetric matrix
    % n is the shape of the symmetric matrix
    % vec has size (N, n(n-1)/2), where N is the number of subjects
    
    N = size(vec, 1);  % Number of subjects
    harmonized_matrix = zeros(N, n, n);  % Initialize the output matrix

    for i = 1:N
        % Create a symmetric matrix for each subject
        tmp = zeros(n);
        inds = triu(true(size(tmp)), 1);  % Indices of upper triangular part (excluding diagonal)
        tmp(inds) = vec(i, :);  % Assign upper triangular part
        tmp = tmp + tmp.';  % Mirror to lower triangular part
        harmonized_matrix(i, :, :) = tmp;
    end
end