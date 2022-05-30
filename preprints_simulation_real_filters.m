%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Prints the results of the files named:
% [method_name]_simulation_real[_white]?.m
% Instructions: load the output file of simulation then run
% this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result_dir = ["./res/"];
simulation_properties = ["filter_1_response_16secs", "filter_2_response_16secs", "filter_3_response_16secs"];
getFullDir = @(file) append(file.folder, '\', file.name);
getInfo = @(error) sprintf("\nmean: %f, std: %f", mean(error), std(error));

estimator_type = ["MLE_WHITE", "BEAMFORMER", "FISHER_SCORING", "FISHER_SCORING_ORIG"];
i=1;
for estim = estimator_type
    est_theta = [];
    est_error_cyclic = [];
    est_error_MSPE = [];
    est_alphas = [];
    re_errors = [];
    re_thetas =[];
    sn = [];
    for prop = simulation_properties
        d = dir(append(result_dir, estim, '_simulation_real_white_results_', prop, '*'));
        if(isempty(d)), fprintf(2, append("property: ", prop, " not found.\n"));end
        load(getFullDir(d(end)));
        est_theta = [est_theta; estimated_theta];
        est_error_cyclic = [est_error_cyclic; estimated_error_cyclic];
        est_error_MSPE = [est_error_MSPE; estimated_error_MSPE];
        est_alphas = [est_alphas; estimated_alphas];
        re_errors = [re_errors; real_errors];
        re_thetas =[re_thetas; real_thetas];
        sn = [sn; snrs];
    end
    [estimated_error_cyclic, indexes] = getMinimumFrom3dArray(est_error_cyclic);
    estimated_theta = getArrayFromIndex(est_theta ,indexes);
    estimated_error_MSPE = getArrayFromIndex(est_error_MSPE ,indexes);
    estimated_alphas = getArrayFromIndex(est_alphas ,indexes);
    real_errors = getArrayFromIndex(re_errors ,indexes);
    real_thetas = getArrayFromIndex(re_thetas ,indexes);
    snrs = getArrayFromIndex(sn ,indexes);
    res = dir(append(result_dir, estim, '_simulation_real_white_results_all_filters_*.mat'));
    save(append(result_dir, estim, "_simulation_real_white_results_all_filters_", string(length(res)+1)),...
        'estimated_theta','real_thetas','estimated_error_cyclic','estimated_error_MSPE','real_errors','estimated_alphas', 'snrs');
end

function [arr, indexes] = getMinimumFrom3dArray(array)
    arr = zeros(1, length(array));
    indexes = zeros(1, length(array));
    j = 1;
    for col=array
        [v, i] = min(col);
        arr(j) = v;
        indexes(j) = i;
        j = j+1;
    end
end

function [arr] = getArrayFromIndex(array, idx)
    arr = zeros(1, length(array));
    i = 1;
    for col=array
        arr(i) = col(idx(i));
        i = i+1;
    end
end




