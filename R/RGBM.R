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

#All steps in the RGBM technique starting from mechanistic network to generation of final inferred GRN
RGBM = function (E      = matrix(rnorm(100),10,10),K = matrix(0,nrow(E),ncol(E)),
                 g_M    = matrix(1,10,10), tfs = paste0("G",c(1:10)), 
                 targets = paste0("G",c(1:10)), lf = 1, 
                 M = 5000, nu = 0.001, s_f = 0.3, no_iterations = 2, mink = 0, 
                 experimentid = 1, outputpath = "DEFAULT",
                 sample_type = "Exp1_",real = 0) 
{
  Ntfs = length(tfs)
  Ntargets = length(targets)
  
  #Add colnames to E matrix and K matrix
  if (is.null(colnames(E)))
  {
    colnames(E) <- paste0("G",c(1:10))
  }
  if (is.null(colnames(K)))
  {
    colnames(K) <- paste0("G",c(1:10))
  }
  
  #Order the mechanistic network in the order of the names of the Tfs and the Targets
  if (is.null(colnames(g_M)))
  {
    colnames(g_M) <- targets;
  }
  else
  {
    g_M <- g_M[,targets];
  }
  if (is.null(rownames(g_M)))
  {
    rownames(g_M) <- tfs;
  }
  else{
    g_M <- g_M[tfs,];
  }
  
  #Get the colids of the tfs active for each target gene from the mechanisitic network
  tf_colids <- foreach(i = 1:Ntargets, .inorder = TRUE, .combine = "cbind") %dopar% {
    result <- rep(0,Ntfs)
    names(result) <- tfs;
    temp <- rownames(g_M)[which(as.numeric(as.vector(g_M[,targets[i]]))>0)];
    result[temp] <- 1
    result
  }
  rm(g_M)
  print("Identified colids for active tfs");
  
  #Get the list of active transcription factors for each target gene
  active_tfs <- c()
  for (i in 1:Ntargets)
  {
    active_tfs <- c(active_tfs,length(which(as.numeric(as.vector(tf_colids[,i]))>0)));
  }
  df_colids <- get_colids(tf_colids,active_tfs,tfs,targets,Ntfs,Ntargets);
  rm(tf_colids)
  print("Identified list of active tfs for each gene")
  
  # First run of GBM + Refinement for the mechanistic set of TFs and Targets
  if (sum(df_colids)==(Ntfs*Ntargets))
  {
    A_first <- first_GBM_step(E,K,tfs,targets,Ntfs,Ntargets,lf,M,nu,s_f,no_iterations)
  }
  else
  {
    A_temp <- matrix(0,nrow=Ntfs,ncol=Ntargets)
    for (iter in 1:no_iterations)
    {
      temp <- second_GBM_step(E,K,df_colids,tfs,targets,Ntfs,Ntargets,lf,M,nu,s_f)
      A_temp <- A_temp + temp
      rm(temp)
    }
    A_temp <- A_temp/no_iterations;
    A_first <- null_model_refinement_step(E,A_temp,K,tfs,targets,Ntfs,Ntargets);
    rm(A_temp)
  }
  gc()
  print("Performed first step of model building")
  
  #Performed the proposed optimal Tf selection using variable importance and detect isolated nodes
  # followed by core GBM model and the refinement step
  if (!dir.exists(outputpath) && outputpath != "DEFAULT")
  {
    dir.create(outputpath)
  }
  if (outputpath == "DEFAULT")
  {
    outputpath = tempdir()
    print(paste("Intermediate GRN and Variable Importance Curve in:", outputpath))
  }
  A <- regularized_GBM_step(E,A_first,K,tfs,targets,Ntfs,Ntargets,lf,M,nu,s_f,experimentid,outputpath,sample_type,mink,real)
  print("Built the final model")
  return(A)
}
