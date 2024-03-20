# QCM-Application-NIC
three data saved in data file folder
step 1: run description0.m to get Tables 1-2; 
step 2: run 2TAR.m and save results in new_new_AUD_USD_cvm_test.mat, new_new_NZD_USD_cvm_test.mat, and new_new_CAD_USD_cvm_test.mat, and Table 4;
step 3: run 3QCM.m and save results in final_CAD_USD_quantiles.mat, final_AUD_USD_quantiles.mat, final_NZD_USD_quantiles.mat,final_CAD_USD_moments.mat, final_AUD_USD_moments.mat, final_NZD_USD_moments.mat;
step 4: run test_moments_skew_kurt_frontier.m to verify the qcm is satisfied skewness-kurtosis frontier;
step 5: run test_moment_by_studentt_test4.m to get Table 3;
step 6: run plot_moment6.R to get Figure 6
step 7: run generate_cv_bandwidth.R and save results in nic_bandwidth.mat;
step 9: run plot_nic_and_bandwidth5.R to get Figure 7 and Table 5;

# QCM-Application-VaR
step 1: run moving_window_qcm.m to get VaR forecast and save results in final_AUD_USD_moments201-500, final_NZD_USD_moments201-500, and final_USD_CAD_moments201-500;
step 2: test_VaR.m to get Table 6;
