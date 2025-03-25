function sc = read_sc(dataframe, C_measure)
    % Number of rows in dataframe
    numRows = height(dataframe);
    
    % Preallocate the sc array
    sc = zeros(numRows, int16(85 * 84 / 2));
    
    % Loop over each row in the dataframe
    for i = 1:numRows
        % Load the .mat file
        matFile = load(dataframe.conn_mat{i});
        
        % Select the matrix based on C_measure
        if strcmp(C_measure, 'C')
            tmp = matFile.C;
        else
            tmp = matFile.aC;
        end
        
        % Vectorize the upper triangular part, excluding the diagonal
        sc(i, :) = tmp(triu(true(size(tmp)), 1));
    end
end