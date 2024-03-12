# QCM
three data saved in data file folder
step 1: run description0.m to get the description results
step 2: run 1exchange_rate_tail_index.m to get the exchange rate series.
step 3: run 2TAR.m and save results in new_new_AUD_USD_cvm_test.mat, new_new_NZD_USD_cvm_test.mat, and new_new_CAD_USD_cvm_test.mat;
step 4: run 3QCM.m and save results in final_CAD_USD_quantiles.mat, final_AUD_USD_quantiles.mat, final_NZD_USD_quantiles.mat;
step 5: run new_QCM.m to get final_CAD_USD_moments.mat, final_AUD_USD_moments.mat, final_NZD_USD_moments.mat;
step 5: run plot_moment6.R to get Figure 6
step 6: run generate_cv_bandwidth.R and save results in nic_bandwidth.mat;
step 7: run plot_nic_and_bandwidth5.R to get Figure 7
step 8: run moving_window_qcm.m to get VaR forecast
step 9: test_VaR.m to get Table 6
