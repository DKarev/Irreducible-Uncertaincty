function [results_V, results_VTU, results_VRU, results_VTURU, results_VTURU_irr] = model_comparison(data)
    
    tbl = data2table(data,1,1); % don't standardize, exclude timeouts
    
    formula = 'C ~ -1 + V + (-1 + V|S)';
    results_V = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + VTU + (-1 + VTU|S)';
    results_VTU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + V + RU + (-1 + V + RU|S)';
    results_VRU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + V + RU + VTU + (-1 + V + RU + VTU|S)';
    results_VTURU = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + V + RU + VTU + RU_irr + (-1 + V + RU + VTU + RU_irr|S)';
    results_VRU_irr = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + V + RU + VTU + VTU_irr + (-1 + V + RU + VTU + VTU_irr|S)';
    results_VTU_irr = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
    
    formula = 'C ~ -1 + V + RU + VTU + RU_irr + VTU_irr + (-1 + V + RU + VTU + RU_irr + VTU_irr|S)';
    results_VTURU_irr = fitglme(tbl,formula,'Distribution','Binomial','Link','Probit','FitMethod','Laplace', 'CovariancePattern','diagonal')
    
    %save results_glme_fig3_nozscore results_V results_VTU results_VRU results_VTURU results_VTURU_irr
    save results_glme_dim_zscore