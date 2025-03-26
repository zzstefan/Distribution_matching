function [harmonized_fea] = harmonize_DM(sc, site_labels, params, weight, reference_site)
    % Get reference site index and parameters
    ref_index = find(strcmp(site_labels, reference_site));
    if strcmp(reference_site, 'OASIS-3')
        ref_param_idx = 1;
        other_sites = {'ADNI2', 'PREVENT-AD'};
    else  % PREVENT-AD
        ref_param_idx = 3;
        other_sites = {'OASIS-3', 'ADNI2'};
    end
    
    % Initialize output
    harmonized_fea = zeros(size(sc));
    
    % Keep reference site unchanged
    for i = 1:size(sc,2)
        harmonized_fea(ref_index,i) = sc(ref_index,i);
        
        % Process other sites
        for j = 1:length(other_sites)
            site = other_sites{j};
            index = find(strcmp(site_labels, site));
            fea = sc(index, i);
            
            % Get correct parameter index for current site
            if strcmp(reference_site, 'OASIS-3')
                site_param_idx = j + 1;
            else
                site_param_idx = j;
            end
            
            % Gamma CDF
            p = gamcdf(fea, params(i, site_param_idx, 1), params(i, site_param_idx, 2));
            
            % Inverse transformation
            p_inv = zeros(size(fea));
            less_than_weight = p <= weight(i, ref_param_idx);
            greater_than_weight = p > weight(i, ref_param_idx);
            p_inv(less_than_weight) = 0;
            p_inv(greater_than_weight) = gaminv((p(greater_than_weight) - weight(i, ref_param_idx)) / ...
                                              (1 - weight(i, ref_param_idx)), ...
                                              params(i, ref_param_idx, 1), ...
                                              params(i, ref_param_idx, 2));
            
            % Assign to harmonized_fea
            harmonized_fea(index,i) = p_inv;
        end
    end
end