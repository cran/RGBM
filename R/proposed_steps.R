#    This is an implementation of RGBM algorithm for Gene Regulatory Network
#    inference from gene/RNA/miRNA expression E, in form of an R package.
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

#Add names to the adjacency matrix
add_names <- function(A,tfs,targets)
{
  rownames(A) <- tfs;
  colnames(A) <- targets;
  return(A);
}

#Get the TF ids 
get_tf_indices <- function(E,tfs,Ntfs){
  
  tf_indices <- c();
  gene_ids <- colnames(E);
  for (i in 1:Ntfs)
  {
    tf_indices <- c(tf_indices,which(gene_ids==tfs[i]));
  }
  rm(E)
  rm(tfs)
  gc()
  return(tf_indices);
}

#Perform the first step of GBM based inference
first_GBM_step <- function(E,K,tfs,targets,Ntfs,Ntargets,lf,M,nu,s_f,no_iterations)
{
  A_first = matrix(0,nrow=Ntfs,ncol=Ntargets);
  tf_indices <- get_tf_indices(E,tfs,Ntfs)
  
  #Construct the averaged RGBM matrix
  for (iter in 1:no_iterations)
  {
    temp <- GBM(E = E, K = K, tfs = tfs, targets=targets, s_s = 1.0, s_f = s_f, lf = lf, M = M, nu = nu, scale = FALSE, center = FALSE)
    A_first <- A_first + temp;
  }
  A_first <- (A_first/no_iterations);
  
  #Perform the null mutant step
  A_first <- null_model_refinement_step(E,A_first,K,tfs,targets,Ntfs,Ntargets);
  
  #Add names to A matrix
  A_first <- add_names(A_first,tfs,targets);
  rm(E)
  rm(K)
  gc()
  return(A_first);
}

#Normalize the adjacency matrix
normalize_matrix_colwise <- function(A,Ntargets)
{
  for (i in 1:Ntargets)
  {
    if (sum(A[,i])>0)
    {
      A[,i] <- A[,i]/sum(A[,i]);
    }
  }
  A[is.nan(A)] <- 0;
  return(A)
}

#Write the first adjacency matrix and get filenames necessary for selection step 
get_filepaths <- function(A_prev,experimentid=1,outputpath,sample_type)
{
  #Write the adjacency matrix after first_GBM_step
  if (!dir.exists(paste0(outputpath,"/Adjacency_Matrix/")))
  {
    dir.create(paste0(outputpath,"/Adjacency_Matrix/"))
  }
  if (!dir.exists(paste0(outputpath,"/Images/")))
  {
    dir.create(paste0(outputpath,"/Images/"))
  }
  write.table(as.data.frame(A_prev),file=paste(outputpath,"/Adjacency_Matrix/",sample_type,"Adjacency_Matrix_",experimentid,".csv",sep=""),col.names = T,row.names = T);
  filepath <- paste(outputpath,"/Adjacency_Matrix/",sample_type,sep="")
  imagepath <- paste(outputpath,"/Images/",sample_type,sep="")
  adjacency_matrix_path <- paste("Adjacency_Matrix_",experimentid,".csv",sep="");
  filepaths <- as.data.frame(cbind(filepath,imagepath,adjacency_matrix_path))
  rm(A_prev)
  gc();
  return(filepaths)
}

#Select ideal number of transcription factors for each target gene
select_ideal_k <- function(experimentid=1,mink=0,filepath,imagepath,adjacency_matrix_path)
{
  
  readfile <- paste(filepath,adjacency_matrix_path,sep="");
  importance_matrix <- read.table(readfile,header=TRUE)
  
  #Keep track of indices while sorting the importance values
  track_matrix <- matrix(0,nrow=dim(importance_matrix)[1],ncol=dim(importance_matrix)[2])
  
  #Make a 4*4 plot to show importance
  pdf(file=paste(imagepath,"PowerGraphImportance_",experimentid,".pdf",sep=""),width = 15, height = 15, pointsize = 20)
  par( mfrow = c( 4, 4 ) )
  N <- ncol(importance_matrix);
  
  for (i in 1:N)
  {
    sorted_ids <- order(importance_matrix[,i],decreasing = TRUE)
    sorted_vector <- importance_matrix[sorted_ids,i]
    track_matrix[,i] <- sorted_ids;
    if (i<=16)
    {
      plot(sorted_vector, type="o", pch=22, lty=2, col ="blue", cex.lab = 1.3, xlab ="Ordered TFs", ylab = "Unnormalized Importance" )
    }
  }
  dev.off()
  
  #Remove targets with nearly uniform distribution of importance weights
  max_importance <- apply(importance_matrix,2,max);
  remove_targets_part1 <- which(max_importance<=.Machine$double.eps);
  if (length(remove_targets_part1)>0){
    max_importance <- max_importance[-remove_targets_part1];
  }
  sorted_importance <- sort(max_importance)
  order_importance <- order(max_importance)
  
  #Since Max_Importance Follows Exponential Distribution Convert to Normal 
  
  #First perform log transformation 
  log_order_importance <- log(sorted_importance)
  #Perform a linear fitting on these transformed values to get a uniform distribution
  #of the fitted values
  x <- c(1:length(log_order_importance))
  y <- as.numeric(as.vector(log_order_importance))
  mod <- lm(y~x)
  fitted_values <- fitted(mod);

  #Transform fitted values into a probability distribution
  fitted_values <- (fitted_values-min(fitted_values))/(max(fitted_values)-min(fitted_values))
    
  #Transform from uniform to normal using the inverse cdf of normal distribution
  normal_values <- qnorm(fitted_values);
  names(normal_values) <- names(sorted_importance);
  
  #Select isolated nodes based on the Gaussian distribution
  consider_normal_values <- normal_values[!is.infinite(normal_values)]
  mean_nv <- mean(consider_normal_values)
  sd_nv <- sd(consider_normal_values)
  cutoff_nv <- mean_nv - 1.5*sd_nv
  remove_targets_part2 <- which(normal_values<cutoff_nv)
  if (length(remove_targets_part2)==0)
  {
    remove_targets_part2 <- NULL ;
  }
  
  #Targets to consider
  targets_to_consider_list <- setdiff(c(1:N),as.numeric(remove_targets_part1));
  targets_to_consider <- setdiff(targets_to_consider_list,order_importance[remove_targets_part2]);
  
  #Calculate ideal k for each target after considering global importance distribution using concept of angles
  Actual_N <- length(targets_to_consider);
  TF <- nrow(importance_matrix);
  ideal_k_vector <- foreach(i = 1:Actual_N, .combine=c , .inorder = TRUE) %dopar%
  {
    target_id <- targets_to_consider[i];
    sorted_tf_impv <- importance_matrix[track_matrix[,target_id],target_id];
    xcord <- c(1:TF);
    ref_v_y <- sorted_tf_impv[length(sorted_tf_impv)]
    ref_v_x <- TF/TF;
    final_vector <- NULL;
    if (mink==0){
      min_id = 1;
    }else{
      min_id=mink;
    }
    for (j in min_id:(TF-2))
    {
      ref_j_y <- sorted_tf_impv[j];
      ref_j_x <- j/TF;
      temp <- c();
      for (k in ((j+1):(TF-1)))
      {
        ref_k_y <- sorted_tf_impv[k];
        ref_k_x <- k/TF;
        a = c(ref_k_x - ref_j_x, ref_k_y - ref_j_y)
        b = c(ref_k_x - ref_v_x, ref_k_y - ref_v_y)
        norm_a <- sqrt(sum(a^2))
        norm_b <- sqrt(sum(b^2))
        numerator <- as.numeric(a%*%b);
        denominator <- (norm_a*norm_b);
        costheta <- numerator/denominator;
        if (costheta<=-0.99999999999999999)
        {
          val <- acos(-1);
        }
        else{
          val <- acos(costheta);
        }
        temp <- c(temp,pi-val);
      }
      select_value <- min(temp);
      select_k<-which.min(temp);
      if ((select_k + j) > (mink + 1))
      {
        final_vector <- rbind(final_vector,cbind(j,(select_k+j),select_value))
      }
      if ((select_k+j)>(mink+1) && select_k==1 && temp==0) {
        break;
      }
    }
    #final_vector <- final_vector[order(final_vector[,2]),];
    min_final_vector_id <- which.min(final_vector[,3]);
    ideal_k <- as.numeric(final_vector[min_final_vector_id,1])-1;
    ideal_k
  }
  rm(importance_matrix)
  gc() 
  #Keep track of the ideal k for each target
  final_k_vector <- rep(0,N);
  for (i in 1:N)
  {
    index <- which(targets_to_consider==i);
    if (length(index)>0)
    {
      final_k_vector[i] <- ideal_k_vector[index];
    }
  }
  return(final_k_vector)
}

#Get the list of transcription factors for each target gene
get_colids <- function(A,ideal_k,tfs,targets,Ntfs,Ntargets)
{
  i = 0
  df_colids = foreach(i = 1:Ntargets, .inorder = TRUE, .combine = "cbind") %dopar% {
    temp <- rep(0,Ntfs);
    names(temp) <- tfs;
    if (ideal_k[i]>0)
    {
      result <- order(as.vector(A[,i]),decreasing=TRUE)
      if (ideal_k[i]<Ntfs)
      {
        result[(ideal_k[i]+1):Ntfs] = 0;
      }
      temp[result] <- 1;
    }
    temp
  }
  rownames(df_colids) <- tfs;
  colnames(df_colids) <- targets;
  rm(A)
  gc()
  return(df_colids);
}

#Construct the second-step GBM
second_GBM_step <- function(E,K,df_colids,tfs,targets,Ntfs,Ntargets,lf,M,nu,s_f)
{
  i = 0
  A_temp <- foreach(i = 1:Ntargets, .inorder = TRUE, .combine = "cbind") %dopar% {
    predictedI  = targets[i];
    predictorsI = tfs[which(as.numeric(df_colids[,predictedI])>0)];
    predictorsI = setdiff(predictorsI,predictedI);
    validE = K[,predictedI] == 0
    result = rep(0,Ntfs);
    names(result) <- tfs;
    len_predictors <- length(predictorsI);
    if (len_predictors>1)
    {
      model = GBM.train(X.train = E[validE,predictorsI],
                         Y.train = as.vector(E[validE,predictedI]),
                         M.train = M,
                         nu = nu,
                         s_s = 1.0,
                         s_f = s_f,
                         lf = lf)
      importance = model$importance;
      importance[is.na(importance)] = 0;
      result[predictorsI] = importance;
    }
    else if (len_predictors==1){
      result[predictorsI] <- 1;
    }
    result
  }
  rm(E)
  rm(K)
  gc()
  A_temp <- add_names(A_temp,tfs,targets);
  return(A_temp)
}

#Apply row deviation as a post-post-processing step
apply_row_deviation <- function(A,Ntfs,Ntargets){
   
   #Apply the row variance 
   s_r  = apply(A,1,sd)
   SC = matrix(rep(s_r,Ntargets),Ntfs,Ntargets)
   A  = A * SC
   return(A);
}

#Get indices for KO experiments
get_ko_experiments <- function(K)
{
  ko.experiments = which(rowSums(K)==1 & apply(K,1,max)==1)
  return(ko.experiments);
}

#Perform null model score if knockout/knockdown information available
z_score_effect <- function(E,K,tfs,targets,Ntfs,Ntargets)
{
  #Build the S2 matrix
  ko.experiments = get_ko_experiments(K)
  S2 = matrix(1,Ntfs,Ntargets)
  rownames(S2) <- tfs;
  colnames(S2) <- targets;
  if (length(ko.experiments)>1) {
    E.ko = E[ko.experiments,]
    rm(E)
    gc()
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
  }
  return(S2);
}

#Apply the null model refinement step
null_model_refinement_step <- function(E,A,K,tfs,targets,Ntfs,Ntargets)
{
  ko.experiments <- get_ko_experiments(K);
  S2 <- z_score_effect(E,K,tfs,targets,Ntfs,Ntargets);
  A <- A*S2;
  rm(E)
  rm(K)
  gc()
  return(A)
}

#Remember the past to help not removing true positives
consider_previous_information <- function(A,A_prev,real)
{
  #Utilize Past Information also to not remove true positives
  A_prev[A_prev==0] <- .Machine$double.eps;
  if (real==0)
  {
    A_prev <- transform_importance_to_weights(A_prev);
  }
  A[A==0] <- .Machine$double.eps;
  if (real==0)
  {
    A <- transform_importance_to_weights(A);
    epsilon <- 1/log(1/.Machine$double.eps);
  }
  else
  {
    epsilon <- .Machine$double.eps;
  }
  A_final <- 2*A*A_prev/(A+A_prev);
  A_final[A_final<=epsilon] <- 0.0;
  A_final[is.nan(A_final)] <- max(A_final[!is.nan(A_final)])
  return(A_final);
}

#Convert the importance scores into edge weights
transform_importance_to_weights <- function(A)
{
  i = 0
  d = ncol(A);
  A_final <- foreach (i = 1:d, .inorder = TRUE, .combine = cbind) %dopar% 
  {
    importance <- as.numeric(as.vector(A[,i]));
    weights <- (1/log(1/importance));
    weights
  }
  return(A_final)
}

#Regulate the size of the regulon by removing out non-significant regulations per TF
regulate_regulon_size <- function(A)
{
  ntfs <- nrow(A);
  for (i in 1:ntfs)
  {
    mean_imp <- mean(A[i,A[i,]>0]);
    sd_imp <- sd(A[i,A[i,]>0]);
    A[i,A[i,]<=(mean_imp+sd_imp)] <- 0;
  }
  return(A)
}

#Perform the proposed regularized GBM
regularized_GBM_step <- function(E,A_prev,K,tfs,targets,Ntfs,Ntargets,lf,M,nu,s_f,experimentid,outputpath,sample_type,mink=0,real=0)
{
  #Get filepaths and allocate them
  filepaths <- get_filepaths(A_prev,experimentid,outputpath,sample_type);
  filepath <- as.character(filepaths$filepath);
  imagepath <- as.character(filepaths$imagepath);
  adjacency_matrix_path <- as.character(filepaths$adjacency_matrix_path);
  
  #Number of columns to select for refinement
  ncol_select <- select_ideal_k(experimentid,mink,filepath,imagepath,adjacency_matrix_path)
  
  #Get the set of trancscription factors active for each target gene
  df_colids <- get_colids(A_prev,ncol_select,tfs,targets,Ntfs,Ntargets);
  
  #Create the refined and regularized GBM model
  A_temp <- second_GBM_step(E,K,df_colids,tfs,targets,Ntfs,Ntargets,lf,M,nu,0.3);
  
  #Apply the null model recitification step
  A <- null_model_refinement_step(E,A_temp,K,tfs,targets,Ntfs,Ntargets);
  rm(A_temp)
  gc()
  
  #Apply the row wise deviation
  A_new <- apply_row_deviation(A,Ntfs,Ntargets);
  rm(A)
  gc()
  
  #Normalize matrix colwise to regularize variations for each target
  A_latest <- normalize_matrix_colwise(A_new,Ntargets);
  rm(A_new)
  gc()
    
  #Remember past information
  A_revised <- consider_previous_information(A_latest,A_prev,real);
  
  #Perform post-processing for real data to regulate size of regulon
  if (real!=0)
  {
    A_revised <- regulate_regulon_size(A_revised);
  }
  
  #Add names to A_final
  A_final <- add_names(A_revised,tfs,targets);
  
  return(A_final)
}


