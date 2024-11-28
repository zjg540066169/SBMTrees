BMTrees_prediction = function(X_train, Y_train, Z_train, subject_id_train, X_test, Z_test, subject_id_test, model = c("BMTrees", "BMTrees_R", "BMTrees_RE", "mixedBART"), binary = FALSE, nburn = 3000L, npost = 4000L, skip = 1L, verbose = TRUE, seed, tol = 1e-20, resample = 5, ntrees = 200, pi_CDP = 0.99){
  n_train = dim(X_train)[1]
  n_test = dim(X_test)[1]
  
  X = rbind(X_train, X_test)
  Y = c(Y_train, Y_test)
  
  
}





#BMTrees_mcmc(NumericMatrix X, NumericVector Y, Nullable<NumericMatrix> Z, CharacterVector subject_id, LogicalVector obs_ind, bool binary = false, long nburn = 0, long npost = 3, bool verbose = true, bool CDP_residual = false, bool CDP_re = false, Nullable<long> seed = R_NilValue, double tol = 1e-40, long ntrees = 200, int resample = 0, double pi_CDP = 0.99)