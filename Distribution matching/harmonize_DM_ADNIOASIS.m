function [harmonized_fea] = harmonize_DM_ADNIOASIS(sc, data, params, weight)
    s1_index = find(strcmp(data.SITE, 'OASIS-3'));
    harmonized_fea = zeros(size(sc));
    for i =1:size(sc,2)
        harmonized_fea(s1_index,i) = sc(s1_index,i);

        % Process other sites
        sites = {'ADNI2'};
        for j = 1:length(sites)
            site = sites{j};
            index = find(strcmp(data.SITE, site));
            fea = sc(index, i);

            % Gamma CDF
            p = gamcdf(fea, params(i, j+1, 1), params(i, j+1, 2)); % MATLAB uses scale as the second parameter

            % Inverse transformation
            p_inv = zeros(size(fea));
            less_than_weight = p <= weight(i, 1);
            greater_than_weight = p > weight(i, 1);
            p_inv(less_than_weight) = 0;
            p_inv(greater_than_weight) = gaminv((p(greater_than_weight) - weight(i, 1)) / (1 - weight(i, 1)), params(i, 1, 1), params(i, 1, 2)); % MATLAB uses scale as the second parameter

            % Assign to harmonized_fea
            harmonized_fea(index,i) = p_inv;
        end
    end

end