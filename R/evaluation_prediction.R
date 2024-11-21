bias = function(Y, Y_O){
  apply(Y, 1, function(x) mean(x - Y_O))
}

absolute_bias = function(Y, Y_O){
  apply(Y, 1, function(x) mean(abs(x - Y_O)))
}

MSE = function(Y, Y_O){
  apply(Y, 1, function(x) mean((x - Y_O)^2))
}

CI_coverge = function(Y, Y_O, symmetry = T){
  mean(sapply(1:length(Y_O), function(i){
    if(symmetry)
      return(quantile(Y[,i], 0.025) <= Y_O[i] & quantile(Y[,i], 0.975) >= Y_O[i])
    else{
      y = Y[,i]
      class(y) = "mcmc"
      ci = coda::HPDinterval(y, prob = 0.95)
      return(ci[1] <= Y_O[i] & ci[2] >= Y_O[i])
    }
  }))
}

Width = function(Y, Y_O, symmetry = T){
  mean(sapply(1:length(Y_O), function(i){
    if(symmetry)
      return(quantile(Y[,i], 0.975) - quantile(Y[,i], 0.025))
    else{
      y = Y[,i]
      class(y) = "mcmc"
      ci = coda::HPDinterval(y, prob = 0.95)
      return(ci[2] - ci[1])
    }
  }))
}


