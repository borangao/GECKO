###R function####
##Data Processing and Jackknife
#Z1<-read.table("summary_-0.15_53_z1.txt",header=T)
#Z2<-read.table("summary_-0.15_53_z2.txt",header=T)
#n1in = 
#n2in =
#nsin = 
#LDscorein
#Weightin = 

#LDscorein<-read.table("gnova_annot_ldscore.csv.gz",header=T)
#Annotation_Lable_test<-read.table("wtccc_ctrl.annot.gz",header=T)


#GECKO(1400,1538,0,as.matrix(unname(LDscore_ordered_L2)),as.matrix(Zscore),as.matrix(unname(Annotation_Lable_test)),T,T,F)
#test_out<-GECKO_R(Z1,Z2,n1in,n2in,nsin,LDscorein,Weightin,Fix_Vein,Test,Annotation_Lable_test)
#' @author Boran Gao, \email{borang@umich.edu}
#' @title
#' GECKO
#' @description
#' Genetic and Environmental Covariance Estimation by Composite-Likelihood Optimization
#'
#' @param Z1 summary statistics for trait 1. The summary statistic file should contain SNP ID, Major allele, and Z statistics with column name of SNP, A1 and Z
#' @param Z2  summary statistics for trait 2. The summary statistic file should contain SNP ID, Major allele, and Z statistics with column name of SNP, A1 and Z
#' @param n1in  Number of samples in study for trait 1
#' @param n2in  Number of samples in study for trait 2
#' @param nsin  Number of  overlapping samples of trait 1 and trait 2. If there is no overlapping sample, then the nsin is zero
#' @param LDscorein LD score for SNPs included in the study. The LD score should be pre-processed by LDSC
#' @param Weightin  Whetehr to include weight for composite-likelihood. Suggested to be true to improve the estimation accuracy
#' @param Test Whethter to test for the genetic and environmental covariance. 
#' @param Annotation_labelin Default is NULL, only need to be included if annotation-stratified genetic and environmental covariance are to be estimated
#' @return List of model parameters including heritability, genetic covariance, environmental covariance, genetic correlation, and environmental correlation estimate and corresponding standard error and P value
#' @export
GECKO_R<-function( Z1,Z2,n1in,n2in,nsin,LDscorein,Weightin,Fix_Vein,Test,Annotation_labelin = NULL){

  Z2_ordered<-Z2[match(Z2$SNP,Z1$SNP),]
  
  Zscorein<-cbind(Z1$Z,Z2_ordered$Z)
  LDscore_ordered<-LDscorein[match(LDscorein$SNP,Z1$SNP),]
  LDscore_ordered_L2<-c()
  Annotation_label<-matrix()
  Annotation_label_name<-c()
  if(missing(Annotation_labelin)){
    Annotation_label<-as.matrix(rep(1,dim(LDscorein)[1]),ncol=1)
    Annotation_label_name<-c("All")
    LDscore_ordered_L2<-as.matrix(LDscore_ordered$L2)
  }else{
    Annotation_label<-as.matrix(Annotation_labelin)
    Annotation_label_name<-colnames(Annotation_labelin)
    LDscore_ordered_L2<-as.matrix(LDscore_ordered[,7:(dim(Annotation_labelin)[2]+6)])
  }

  
  
  if(Test == T){
    result = GECKO(n1in,n2in, nsin,LDscore_ordered_L2,Zscorein,Annotation_label,Weightin,Fix_Vein,FALSE)
    Est<-c()
    for(j in 1:(dim(Annotation_label)[2]+1)){
      
      j_th_res<-result[,,j][upper.tri(result[,,j],diag=T)]
      j_th_res<-c(j_th_res,j_th_res[2]/sqrt(j_th_res[1]*j_th_res[3]))
      Est<-c(Est,j_th_res)
    }
    
    
    
    label_index<-sort(rep(1:100,ceiling(dim(Zscorein)[1]/100))[1:dim(Zscorein)[1]])
   # tic<-proc.time()
    Jack_result<-list()
    for(i in 1:100){
      Zscorein_subset<-as.matrix(Zscorein[label_index!=i,])
      Annotation_labelin_subset<-as.matrix(Annotation_label[label_index!=i,])
      LDscorein_subset<-as.matrix(LDscore_ordered_L2[label_index!=i,])
      result_subset = GECKO(n1in,n2in, nsin,LDscorein_subset,Zscorein_subset,Annotation_labelin_subset,Weightin,Fix_Vein,Test)
      Jack_result[[i]] = result_subset
      cat(i)
    }
   # toc<-proc.time()-tic
   # print(toc)
      Jack_result_summary<-c()
      for(i in 1:100){
      i_th_iter<-c()
      for(j in 1:(dim(Annotation_label)[2]+1)){
        
        j_th_res<-Jack_result[[i]][,,j][upper.tri(Jack_result[[i]][,,j],diag=T)]
        j_th_res<-c(j_th_res,j_th_res[2]/sqrt(j_th_res[1]*j_th_res[3]))
        
        i_th_iter<-c(i_th_iter,j_th_res)
      }
      Jack_result_summary<-rbind(Jack_result_summary,i_th_iter)
      }
    Jackknife_sd<-sqrt(apply(Jack_result_summary,2,var)*99*99/100)
    Jackknife_pvalue<-2*pnorm(-abs(Est/Jackknife_sd))  
    
    output<-rbind(Est,Jackknife_sd,Jackknife_pvalue)
    dimnames(output)[[1]]<-c("Estimates","Standard_Error","P_value")
    
    
    if(length(Annotation_label_name)==1){
      dimnames(output)[[2]]<-c("Env_Var_Comp1","Env_Covar","Env_Var_Comp2","Env_Correlation","Gen_Var_Comp1","Gen_Covar","Gen_Var_Comp2","Gen_Corr")
    }else{
      annotataion_col_names<-c()
      for(j in 1:length(Annotation_label_name)){
        annotataion_col_names<-c(annotataion_col_names,paste0(Annotation_label_name[j],c("_Gen_Var_Comp1","_Gen_Covar","_Gen_Var_Comp2","_Gen_Corr")))
      }
      dimnames(output)[[2]]<-c("Env_Var_Comp1","Env_Covar","Env_Var_Comp2","Env_Correlation",annotataion_col_names)
      
    }
       return(output)    
  }else{
    result = GECKO(n1in,n2in, nsin,LDscore_ordered_L2,Zscorein,Annotation_label,Weightin,Fix_Vein,Test)
    Est<-c()
    for(j in 1:(dim(Annotation_label)[2]+1)){
      
      j_th_res<-result[,,j][upper.tri(result[,,j],diag=T)]
      j_th_res<-c(j_th_res,j_th_res[2]/sqrt(j_th_res[1]*j_th_res[3]))
      Est<-c(Est,j_th_res)
    }
    Est<-matrix(Est,nrow = 1)
    
    dimnames(Est)[[1]]<-list("Estimates")
    
    if(length(Annotation_label_name)==1){
      dimnames(Est)[[2]]<-c("Env_Var_Comp1","Env_Covar","Env_Var_Comp2","Env_Correlation","Gen_Var_Comp1","Gen_Covar","Gen_Var_Comp2","Gen_Corr")
    }else{
      annotataion_col_names<-c()
      for(j in 1:length(Annotation_label_name)){
        annotataion_col_names<-c(annotataion_col_names,paste0(Annotation_label_name[j],c("_Gen_Var_Comp1","_Gen_Covar","_Gen_Var_Comp2","_Gen_Corr")))
      }
      dimnames(Est)[[2]]<-c("Env_Var_Comp1","Env_Covar","Env_Var_Comp2","Env_Correlation",annotataion_col_names)
      
    }
        return(Est)
  }
}

#write.table(test_out,"test.txt",col.names = T,row.names = F,quote=F)


