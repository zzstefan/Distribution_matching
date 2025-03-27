% Adjust the code, without the large datatable, we now use directly use the
% site_labels as input;

function [weight, params] = fit_dist_new(conn_mat, site_labels, fea_index)
    % Nested function for mixed distribution
    function p = mixed_distribution(x, amplitude, shape, scale)
        p = amplitude * (x == 0) + (1 - amplitude) * gampdf(x, shape, scale);
    end

    %SITE = {'Site 1', 'Site 2', 'Site 3'};
    SITE = {'OASIS-3', 'ADNI2', 'PREVENT-AD'};
    weight = zeros(1, length(SITE));
    params = zeros(length(SITE), 2);
    
    % Create figure and subplots
    %fig = figure;
    % fig = figure('Visible', 'off');
    % set(fig, 'Position', [0, 0, 450, 150]);
    for i = 1:length(SITE)
        site = SITE{i};
        index = find(strcmp(site_labels, site));
        sc = conn_mat(index, fea_index);
        num_zero = sum(sc == 0);
        weight(i) = num_zero / length(sc);

        % x_values = linspace(min(conn_mat, [], 'all'), max(conn_mat, [], 'all'), 1000);
        % Plot the data histogram
        % subplot(1, 3, i);
        % histogram(sc, 40, 'Normalization', 'pdf');
        % hold on;
        
        % Plot the fitted mixed distribution
        num_nonzero = sum(sc ~= 0);
        amplitude = num_zero / length(sc);
        if num_nonzero >=2
            gamma_params = gamfit(sc(sc ~= 0), [], []); % floc equivalent in MATLAB
            params(i, :) = [gamma_params(1), gamma_params(2)];
            % plot(x_values, mixed_distribution(x_values, amplitude, gamma_params(1), gamma_params(2)), 'r-');
        elseif num_nonzero == 1
            % MATLAB doesn't have direct equivalent of f0, need to handle this case manually
            % Example: fixed shape parameter
            shapeParam = 22.42; % Example value
            params(i, 1) = shapeParam;
            scaleParam = mean(sc(sc ~= 0)) / shapeParam; % Example estimation for scale
            params(i, 2) = scaleParam;
            % plot(x_values, mixed_distribution(x_values, amplitude, shapeParam,scaleParam), 'r-');
        else
            params(i, :) = [NaN, NaN];
        end
        
        % title(sprintf('Site %d', i));
        % hold off;
    end
    % Adjust layout and save the figure
    %saveas(fig, sprintf('./Result/fitted_distribution_all_matlab/fea_%d.png', fea_index), 'png');
    % exportgraphics(fig, sprintf('./Result/fitted_distribution_all_matlab/fea_%d.png', fea_index), 'Resolution', 300);
    % close(fig);
end