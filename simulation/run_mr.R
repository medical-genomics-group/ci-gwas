
library("igraph");library("bnlearn")
library("MRPRESSO");library("mixtools")
library("optparse");library("limma")
library("MendelianRandomization")
library("Matching");library("rmatio")
library("Zelig");library("Metrics")
library("CBPS");library("pcalg")
library("MatchIt");library("cobalt")
library("fastmatch");library("stdReg")
library("Matrix");library("readr")
library("dplyr");library("cause")
library("parallel");library("ParallelPC")
library("ff", quietly = TRUE)
library("ivreg");library("mclust")
library("tictoc")

source("./MRfunctions.R")

myargs = commandArgs(trailingOnly=TRUE)
alpha_e = as.numeric(myargs[1])
id = as.numeric(myargs[2])

set.seed(id)
n=16000
SNP=1600


load(file=paste( "./sim_small_effects/test_n", toString(n), "_SNP_", toString(SNP),"_Rep_",toString(id) , ".RData", sep = "" ))

Tr=10
#y_pa=5
nL=2
deg=3
prob_pleio=0.2
pq=SNP+Tr+nL

Bg=myDAG3$G[(SNP+nL+1):pq,(SNP+nL+1):pq]
Bg=t(Bg)
Bg = igraph::graph_from_adjacency_matrix(t(abs(Bg)>0))
#plot(igraph::simplify(Bg), vertex.color="green")
B_distances = igraph_directed_distances(Bg)
#print(B_distances)

Y=myDAG3$x[,(SNP+nL+1):pq]
###################################################################################
# Create the input for MR
all_mr_res = list()
df = dat
ivs = colnames(df)[grepl("X",colnames(df))]
phenos = colnames(df)[grepl("Y",colnames(df))]
num_ivs = length(ivs)
p = length(phenos)
# Get all IV-phenotype associations
GWAS_Ps = matrix(1,num_ivs,p, dimnames = list(ivs,phenos))
GWAS_effects = matrix(0,num_ivs, p, dimnames = list(ivs,phenos))
GWAS_ses = matrix(0,num_ivs, p, dimnames = list(ivs,phenos))
GWAS_Zs = matrix(0,num_ivs, p, dimnames = list(ivs,phenos))

for (pheno in phenos) {
    print(pheno)
    gwas_res = sapply(ivs,run_lm,x=pheno,z=NULL,df = df)
    GWAS_Ps[,pheno] = gwas_res[4,]
    GWAS_effects[,pheno] = gwas_res[1,]
    GWAS_ses[,pheno] = gwas_res[2,]
    GWAS_Zs[,pheno] = gwas_res[3,]
}

print(alpha_e)
G_it = GWAS_Ps < alpha_e
res_path="./mr_res_git_thr"
MMfile_git <- file.path(res_path, paste("git_n", toString(n), "_SNP_", toString(SNP),  "_alpha_", toString(alpha_e), "_it_", toString(id), ".csv", sep = ""))
print(MMfile_git)
   print(G_it)
  write.csv(G_it, file=MMfile_git,row.names = FALSE)
  # Run MR
  # pleio size is set to 100 - no filtering of variants (a standard MR analysis)
  mr_anal_res = list(
    "Egger" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
                                       pleio_size=100,pruned_lists=NULL,func=mr_egger,robust=T),
    "IVW" = run_pairwise_mr_analyses(G_it,GWAS_effects,GWAS_ses,
                                     pleio_size=100,pruned_lists=NULL,func=mr_ivw,robust=T)
  )
  print(mr_anal_res)


  all_mr_res[["egger"]] = rbind(all_mr_res[["egger"]],mr_anal_res$Egger)
  all_mr_res[["ivw"]] = rbind(all_mr_res[["ivw"]],mr_anal_res$IVW)
  
MMfile1 <- file.path(res_path, paste("mr_skeleton_egger_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = ""))##
MMfile2 <- file.path(res_path, paste("mr_skeleton_ivw_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = ""))##
write.csv(all_mr_res[1], file = MMfile1, row.names = FALSE)
write.csv(all_mr_res[2], file = MMfile2, row.names = FALSE)
  # Add MRPRESSO
  mrpresso_res = c()
  #cgauge_mrpresso_res = c()
  try({
    for(tr1 in phenos){
      currivs = rownames(GWAS_Ps)[G_it[,tr1]]
      for(tr2 in phenos){
        if(tr1==tr2){next}
        X = data.frame(E1b=GWAS_effects[currivs,tr1],O1b=GWAS_effects[currivs,tr2],
                       E1sd=GWAS_ses[currivs,tr1],O1sd=GWAS_ses[currivs,tr2])
        try({
          res = mr_presso(BetaOutcome = "O1b", BetaExposure = "E1b", 
                          SdOutcome = "O1sd", SdExposure = "E1sd",data=X,
                          OUTLIERtest=T,
                          DISTORTIONtest = T,
                          NbDistribution = 1000,SignifThreshold = 0.1)
          if(is.na(res$`Main MR results`[2,"P-value"])){
            v = c(tr1,tr2,unlist(res$`Main MR results`[1,]))
          }
          else{
            v = c(tr1,tr2,unlist(res$`Main MR results`[2,]))
          }
          v["GlobalTestP"] = res$`MR-PRESSO results`$`Global Test`$Pvalue
          mrpresso_res = rbind(mrpresso_res,v)
        })
      }
    }
    if(!is.null(dim(mrpresso_res))){
      mrpresso_res = as.data.frame(mrpresso_res)
      for(j in 3:ncol(mrpresso_res)){
        mrpresso_res[[j]] = as.numeric(as.character(mrpresso_res[[j]]))
      }  
    }
  })
  print( mrpresso_res )

  all_mr_res[["mrpresso"]] = rbind(all_mr_res[["mrpresso"]],mrpresso_res)
  MMfile3 <- file.path(res_path, paste("mr_skeleton_mrpresso_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = ""))##
  write.csv(all_mr_res[3], file = MMfile3, row.names = FALSE)


  # Add CAUSE estimates - 
   cause_res = mclapply(phenos,run_cause_on_tr,
                        phenos=phenos,G_it=G_it,GWAS_effects=GWAS_effects,GWAS_ses=GWAS_ses,
                        B_distances = B_distances,
                        mc.cores = 4)
   cause_res_all = c()
   for(tr in cause_res){
     if(length(tr)==1){next}
     cause_res_all= rbind(cause_res_all,tr)
   }
   cause_res_all_df = data.frame(cause_res_all,stringsAsFactors = F)
  for(j in 3:ncol(cause_res_all_df)){
     cause_res_all_df[[j]] = as.numeric(as.character(cause_res_all_df[[j]]))
}
 


  
  
  
all_mr_res[["CAUSE"]] = rbind(all_mr_res[["CAUSE"]],cause_res_all_df)

save(all_mr_res, file=paste( "./mr_res_git_thr/all_mr_res_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_Rep_",toString(id) , ".RData", sep = "")) 

res_path="./mr_res_git_thr"
MMfile4 <- file.path(res_path, paste("mr_skeleton_cause_n", toString(n), "_SNP_", toString(SNP),"_alpha_", toString(alpha_e),"_it_",toString(id) , ".csv", sep = ""))##
#merged_all_mr <- bind_rows(all_mr_res [1], all_mr_res[2],all_mr_res[3],all_mr_res[4])

write.csv(all_mr_res[4], file = MMfile4, row.names = FALSE)
