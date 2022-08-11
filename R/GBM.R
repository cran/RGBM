#    This is an implementation of RGBM algorithm for Gene Regulatory Network
#    inference from gene/RNA/miRNA expression data, in form of an R package.
#    Copyright (C) 2016  Raghvendra Mall

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program, see LICENSE.

GBM = function (E      = matrix(rnorm(100),10,10),
                  K      = matrix(0,nrow(E),ncol(E)),
                  tfs     = paste0("G",c(1:10)),
                  targets = paste0("G",c(1:10)),
                  s_s    = 1,
                  s_f    = 0.3,
                  lf     = 1,
                  M      = 5000,
                  nu     = 0.001,
                  scale  = TRUE,
                  center = TRUE,
                  optimization.stage  = 2) {
  
  S = nrow(E)
  E = scale(E,scale=scale,center=center)
  #Add colnames to E matrix and K matrix
  if (is.null(colnames(E)))
  {
    colnames(E) <- paste0("G",c(1:10))
  }
  if (is.null(colnames(K)))
  {
    colnames(K) <- paste0("G",c(1:10))
  }
  
  Ntfs = length(tfs)
  Ntargets = length(targets)
  V = matrix(0, Ntfs, Ntargets)
  i = 0
  V = foreach(i = 1:Ntargets, .inorder = TRUE, .combine = "cbind") %dopar% {
    predictedI  = targets[i];
    predictorsI = setdiff(tfs, predictedI)
    validE      = K[,predictedI] == 0
    
    model = GBM.train(X.train = E[validE,predictorsI],
                        Y.train = E[validE,predictedI],
                        M.train = M,
                        nu = nu,
                        s_s = s_s,
                        s_f = s_f,
                        lf = lf)
    result = rep(0, Ntfs)
    names(result) = tfs;
    result[predictorsI] = model$importance
    result
  }
  rownames(V) <- tfs;
  colnames(V) <- targets;
  
  if (optimization.stage > 0) {
    # first stage of re-evaluation
    s  = apply(V,1,var)
    S1 = matrix(rep(s,Ntargets),Ntfs,Ntargets)
    V  = V * S1
  }
  
  if (optimization.stage > 1) {
    # second stage of re-evaluation
    ko.experiments = which(rowSums(K)==1 & apply(K,1,max)==1)
    if (length(ko.experiments)>1) {
      S2 = matrix(1,Ntfs,Ntargets)
      rownames(S2) <- tfs;
      colnames(S2) <- targets;
      E.ko = E[ko.experiments,targets]
      for (tf in tfs) {
        for (target in targets) {
          if (sum(K[ko.experiments,tf]==1)>0)
          {
            avg.ko   = mean(E.ko[K[ko.experiments,tf]==1,target])
            avg.n.ko = mean(E.ko[K[ko.experiments,tf]==0,target])
            std.dev  = sd(E.ko[,target])
            if (std.dev>0) {
              S2[tf,target]  = abs(avg.ko-avg.n.ko)/std.dev
            }
          }
        }
      }
      V = V * S2
    }}
  
  # prepare the final adjacency matrix
  colnames(V) = targets;
  rownames(V) = tfs;
  return(V)
}
