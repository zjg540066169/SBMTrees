library(tidyverse)

bias = function(X, X_O, R, overall = F, posterior = F){
  if(any(class(X_O) == "data.frame")){
    X_O = as.matrix(X_O)
  }
  if(posterior == F){
    if(all(class(X) == "numeric")){
      return(mean((X - X_O)[R]))
    }else{
      if(overall){
        return(mean((X - X_O)[R]))
      }else{
        sapply(1:dim(X)[2], function(i){
          mean((X[,i] - X_O[,i])[R[,i]])
        })
      }
    }
  }else{
    if(length(dim(X)) == 2){
      return(mean((colMeans(X) - X_O)[R]))
    }else{
      if(overall){
        return(mean((apply(X, c(2,3), median) - X_O)[R]))
      }else{
        X = apply(X, c(2,3), median)
        sapply(1:dim(X)[2], function(i){
          mean((X[,i] - X_O[,i])[R[,i]])
        })
      }
    }
  }
}


mse = function(X, X_O, R, overall = F, posterior = F){
  if(any(class(X_O) == "data.frame")){
    X_O = as.matrix(X_O)
  }
  if(posterior == F){
    if(all(class(X) == "numeric")){
      return(mean(((X - X_O)^2)[R]))
    }else{
      if(overall){
        return(mean(((X - X_O)^2)[R]))
      }else{
        sapply(1:dim(X)[2], function(i){
          mean(((X[,i] - X_O[,i])^2)[R[,i]])
        })
      }
    }
  }else{
    if(length(dim(X)) == 2){
      return(mean(((colMeans(X) - X_O)^2)[R]))
    }else{
      if(overall){
        #return(mean(((apply(X, c(2,3), mean) - X_O)^2)[R]))
        return(mean(apply(X, 1, function(x) mean(((x - X_O)[R])^2))))
      }else{
        # X = apply(X, c(2,3), mean)
        # sapply(1:dim(X)[2], function(i){
        #   mean(((X[,i] - X_O[,i])^2)[R[,i]])
        # })
        colMeans(sapply(1:dim(X)[3], function(i){
          apply(X[,,i], 1, function(x){
            mean((x - X_O[,i])^2)
          })
        }))
      }
    }
  }
}

mse_summary = function(X, X_O, R, overall = F, posterior = F){
  if(any(class(X_O) == "data.frame")){
    X_O = as.matrix(X_O)
  }
  if(posterior == F){
    if(all(class(X) == "numeric")){
      return(mean(((X - X_O)^2)[R]))
    }else{
      if(overall){
        return(mean(((X - X_O)^2)[R]))
      }else{
        sapply(1:dim(X)[2], function(i){
          mean(((X[,i] - X_O[,i])^2)[R[,i]])
        })
      }
    }
  }else{
    if(length(dim(X)) == 2){
      return(mean(((colMeans(X) - X_O)^2)[R]))
    }else{
      if(overall){
        return(mean(((apply(X, c(2,3), median) - X_O)^2)[R]))
      }else{
        X = apply(X, c(2,3), median)
        sapply(1:dim(X)[2], function(i){
          mean(((X[,i] - X_O[,i])^2)[R[,i]])
        })
      }
    }
  }
}

mse_mean = function(X, X_O, R, overall = F, posterior = F){
  if(any(class(X_O) == "data.frame")){
    X_O = as.matrix(X_O)
  }
  if(posterior == F){
    if(all(class(X) == "numeric")){
      return(mean(((X - X_O)^2)[R]))
    }else{
      if(overall){
        return(mean(((X - X_O)^2)[R]))
      }else{
        sapply(1:dim(X)[2], function(i){
          mean(((X[,i] - X_O[,i])^2)[R[,i]])
        })
      }
    }
  }else{
    if(length(dim(X)) == 2){
      return(mean(apply(X, 1, function(x) mean(((x - X_O)^2)[R]))))
    }else{
      if(overall){
        return(mean(((apply(X, c(2,3), mean) - X_O)^2)[R]))
      }else{
        return(
          return(sapply(1:dim(X)[3], function(i){
            mse_mean(X[,,i], X_O[,i], R[,i], overall = F, posterior = T)
          }))
        )
      }
    }
  }
}



binary_bias = function(X, X_O, R, overall = F, posterior = F){
  if(any(class(X_O) == "data.frame")){
    X_O = as.matrix(X_O)
  }
  if(posterior == F){
    if(all(class(X) == "numeric")){
      return(mean((X != X_O)[R]))
    }else{
      if(overall){
        return(mean((X != X_O)[R]))
      }else{
        sapply(1:dim(X)[2], function(i){
          mean((X[,i] != X_O[,i])[R[,i]])
        })
      }
    }
  }else{
    if(length(dim(X)) == 2){
      return(mean(apply(X, 1, function(x) mean((x != X_O)[R]))))
    }else{
      if(overall){
        return(mean(apply(X, 1, function(x) mean((x != X_O)[R]))))
      }else{
        sapply(1:dim(X)[3], function(i){
          mean(apply(X[,,i], 1, function(x) mean((x[R[,i]] != X_O[R[,i], i]))))
        })
      }
    }
  }
}


binary_criteria = function(X, X_O, R, overall = F, posterior = F){
  if(any(class(X_O) == "data.frame")){
    X_O = as.matrix(X_O)
  }
  if(posterior == F){
    if(all(class(X) == "numeric")){
      return(mean((X != X_O)[R]))
    }else{
      if(overall){
        return(mean((X != X_O)[R]))
      }else{
        sapply(1:dim(X)[2], function(i){
          mean((X[,i] != X_O[,i])[R[,i]])
        })
      }
    }
  }else{
    if(length(dim(X)) == 2){
      return(mean((apply(X, 2, function(y) DescTools::Mode(y)[1]) != X_O)[R]))
    }else{
      if(overall){
        return(mean((apply(X, c(2,3), function(y) DescTools::Mode(y)[1]) != X_O)[R]))
      }else{
        X = apply(X, c(2,3),function(y) DescTools::Mode(y)[1])
        sapply(1:dim(X)[2], function(i){
          mean((X[,i] != X_O[,i])[R[,i]])
        })
      }
    }
  }
}
