%%%
% For this function, we offer two choices
% First is for just unharmonized data or harmonized data
% Second is for combining both unharmonized and harmonized data together
% for comparison;
% single =1, means one column either unharmonized or harmonized
% single = 0, combining them together
%%%
function fit_dist_histogram_ADNIOASIS(conn_mat,har_fea, data, fea_index,single, sex_group)
    % Nested function for mixed distribution
    function p = mixed_distribution(x, amplitude, shape, scale)
        p = amplitude * (x == 0) + (1 - amplitude) * gampdf(x, shape, scale);
    end
    SITE = {'OASIS-3', 'ADNI2'};
    n_site = length(SITE);

    result_folder = sprintf('./fitted_histogram_%s/', sex_group);
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
    end
    

    if single
        
        weight = zeros(1, n_site);
        params = zeros(n_site, 2);

        % Create figure and subplots
        %fig = figure;
        fig = figure('Visible', 'off');
        set(fig, 'Position', [0, 0, 160, n_site*150+20]);
        for i = 1:n_site
            site = SITE{i};
            index = find(strcmp(data.SITE, site));
            sc = conn_mat(index, fea_index);
            num_zero = sum(sc == 0);
            

            x_values = linspace(min(conn_mat, [], 'all'), max(conn_mat, [], 'all'), 1000);
            % Plot the data histogram
            subplot(n_site, 1, i);
            histogram(sc, 40, 'Normalization', 'pdf');
            hold on;

            % Plot the fitted mixed distribution
            num_nonzero = sum(sc ~= 0);
            amplitude = num_zero / length(sc);
            if num_nonzero >= 2
                gamma_params = gamfit(sc(sc ~= 0), [], []); % floc equivalent in MATLAB
                params(i, :) = [gamma_params(1), gamma_params(2)];
                plot(x_values, mixed_distribution(x_values, amplitude, gamma_params(1), gamma_params(2)), 'r-');
            elseif num_nonzero == 1
                % MATLAB doesn't have direct equivalent of f0, need to handle this case manually
                % Example: fixed shape parameter
                shapeParam = 22.42; % Example value
                params(i, 1) = shapeParam;
                scaleParam = mean(sc(sc ~= 0)) / shapeParam; % Example estimation for scale
                params(i, 2) = scaleParam;
                plot(x_values, mixed_distribution(x_values, amplitude, shapeParam,scaleParam), 'r-');
            else
                params(i, :) = [NaN, NaN];
            end

            title(sprintf('Site %d', i));
            % hold off;
        end
        % Adjust layout and save the figure
        %saveas(fig, sprintf('./Result/fitted_distribution_all_matlab/fea_%d.png', fea_index), 'png');
   

        exportgraphics(fig, sprintf(strcat(result_folder,'/DM_fea_%d.png'), fea_index), 'Resolution', 300);
        close(fig);
    else
        fig = figure('Visible', 'off');
        set(fig, 'Position', [0, 0, 320, n_site*150+20]);
        con_data(:,:,1) = conn_mat;
        con_data(:,:,2) = har_fea;
        x_values = linspace(min(con_data(:,fea_index,:),[], "all"), max(con_data(:,fea_index,:),[], "all"), 1000);
        xLimits = [min(con_data(:,fea_index,:),[], "all") max(con_data(:,fea_index,:),[], "all")];
        for i = 1:length(SITE)
            site = SITE{i};
            index = find(strcmp(data.SITE, site));
            for j = 1:2
    
                % Plot the data histogram
                
                sc = con_data(index, fea_index,j);
                if j==1
                    subplot(n_site, 2, 2*i-1);
                else
                    subplot(n_site,2,2*i);
                end
                
                num_zero = sum(sc == 0);
                
                histogram(sc, 40, 'Normalization', 'pdf');
                xlim(xLimits);
                hold on;

                % Plot the fitted mixed distribution
                num_nonzero = sum(sc ~= 0);
                amplitude = num_zero / length(sc);
                if num_nonzero >= 2
                    gamma_params = gamfit(sc(sc ~= 0), [], []); % floc equivalent in MATLAB
                    params(i, :) = [gamma_params(1), gamma_params(2)];
                    plot(x_values, mixed_distribution(x_values, amplitude, gamma_params(1), gamma_params(2)), 'r-');
                elseif num_nonzero == 1
                    % MATLAB doesn't have direct equivalent of f0, need to handle this case manually
                    % Example: fixed shape parameter
                    shapeParam = 22.42; % Example value
                    params(i, 1) = shapeParam;
                    scaleParam = mean(sc(sc ~= 0)) / shapeParam; % Example estimation for scale
                    params(i, 2) = scaleParam;
                    plot(x_values, mixed_distribution(x_values, amplitude, shapeParam,scaleParam), 'r-');
                else
                    params(i, :) = [NaN, NaN];
                end
                if i~=1
                    if j==1
                        title(sprintf('UnHar ADNI2'));
                    else
                        title(sprintf('Har ADNI2'));
                    end
                else
                    title(sprintf('OASIS-3'));
                end
            % hold off;
            end
        end
        % Adjust layout and save the figure
        %saveas(fig, sprintf('./Result/fitted_distribution_all_matlab/fea_%d.png', fea_index), 'png');

        exportgraphics(fig, sprintf(strcat(result_folder,'/DM_fea_%d.png'), fea_index), 'Resolution', 300);
        close(fig);

    end
end