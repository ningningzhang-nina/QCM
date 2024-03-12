# QCM
step 1: generate simulated data following five DGPs:
        run new_DGP1_garch.m and save results in new_DGP1_garch.mat;
        run new_DGP2_studentt.m and save results in new_DGP2_studentt.mat;
        run new_DGP3_garch_snp.m and save results in new_DGP3_garch_snp.mat;
        run new_DGP4_gir_skew_studentt.m and save results in new_DGP4_gir_skew_studentt.mat;
        run DGP5_arma_mx_garch.m and save results in DGP5_arma_mx_garch.mat
step 2: run r_case5_garchsk_estimation.R to get the garchsk model results.
step 3: run plot_simulation_boxplot.R to get Figures 1-5;
step 4: run plot_threshold_boxplot.R to get Figures B1-5;
step 5: run plot_5dgp_means.R to get Figure B6; run plot_5dgp_sd.R to get Figure B7;
