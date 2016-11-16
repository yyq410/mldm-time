testOtherMethods <- function(X, M, subdir, test_number) {
  library("psych")
  library("huge")
  library("ccrepe")
  library("SpiecEasi")

  data <- X
  datar <- X / rowSums(X)
  Mscale <- scale(M, center = TRUE, scale = TRUE)

############################################
# PCC
############################################
pcc_result <- corr.test(as.data.frame(datar))
pcc_pvalue <- pcc_result$p
pcc_corr <- pcc_result$r

pcc_env_result <- corr.test(x=data.frame(datar), y=data.frame(Mscale))
pcc_env_pvalue <- pcc_env_result$p
pcc_env_corr <- pcc_env_result$r

# record all pcc results
pcc_record <- list(pcc_pvalue, pcc_env_pvalue, pcc_corr, pcc_env_corr)

pcc_file <- paste(subdir,'pcc-',as.character(test_number),'.RData',sep="")
save(pcc_record, file=pcc_file)

rm(pcc_result,pcc_pvalue,pcc_corr,pcc_env_result,pcc_env_pvalue,pcc_env_corr)
rm(pcc_record)
print('after pcc:')
gc()
############################################
# Spearman
############################################
spearman_result <- corr.test(as.data.frame(datar), method="spearman")
spearman_pvalue <- spearman_result$p
spearman_corr <- spearman_result$r

spearman_env_result <- corr.test(x=data.frame(datar), y=data.frame(Mscale), method="spearman")
spearman_env_pvalue <- spearman_env_result$p
spearman_env_corr <- spearman_env_result$r

spearman_env <- spearman_env_corr

# record all pcc results
spearman_record <- list(spearman_pvalue, spearman_env_pvalue, spearman_corr, spearman_env_corr)

#print("spearman otu-otu:")
#print(spearman_corr)
#print("spearman ef-otu:")
#print(spearman_env_corr)

scc_file <- paste(subdir,'scc-',as.character(test_number),'.RData',sep="")
save(spearman_record, file=scc_file)

rm(spearman_result,spearman_pvalue,spearman_corr,spearman_env_corr,spearman_env_pvalue,spearman_env_result)
rm(spearman_record)
print('after spearman:')
gc()

############################################
# graphical lasso
############################################
gl <- huge(cbind(Mscale, datar), method="glasso", cov.output=TRUE, nlambda=30)

# record all glasso results
gl_record <- list(gl)

############################################
# mLasso
############################################
#ml <- huge(cbind(Mscale, datar), method="mb", sym="or", nlambda=30)

# record all mLasso results
ml_record <- list()
gl_file <- paste(subdir,'gl-',as.character(test_number),'.RData',sep="")
save(gl_record, file=gl_file)

rm(gl_record, gl)
print('after glasso:')
gc()

############################################
# sparCC
############################################
sparcc_result <- sparcc(data, iter=50, inner_iter=25, th=0.1)
sparcc_cov <- sparcc_result[[1]]

#print("SparCC result:")
#print(sparcc_cov)

# record all sparCC results
sparcc_record <- list(sparcc_cov, sparcc_result)

sparcc_file <- paste(subdir,'sparcc-',as.character(test_number),'.RData',sep="")
save(sparcc_record, file=sparcc_file)

rm(sparcc_record, sparcc_result, sparcc_cov)
print('after sparcc:')
gc()

############################################
# ccrepe
############################################
ccrepe_result <- ccrepe(x=datar, sim.score=cor, sim.score.args=list(method='pearson'))
ccrepe_sim <- ccrepe_result$sim.score
ccrepe_p <- ccrepe_result$p.values

#print("ccrepe result:")
#print(ccrepe_sim)

# record all ccrepe results
ccrepe_record <- list(ccrepe_p, ccrepe_sim)

ccrepe_file <- paste(subdir,'ccrepe-',as.character(test_number),'.RData',sep="")
save(ccrepe_record, file=ccrepe_file)

rm(ccrepe_result, ccrepe_sim, ccrepe_p, datar)
rm(ccrepe_record)
print('after ccrepe:')
gc()

############################################
# spiec-easi glasso
############################################
spiec_glasso <- spiec.easi(data, method='glasso', sel.criterion='stars', nlambda=30)

############################################
# spiec-easi mb
############################################
spiec_mb <- spiec.easi(data, method='mb', sel.criterion='stars', nlambda=30)

# record all spiec-easi results
spiec_record <- list(spiec_glasso, spiec_mb)

spiec_file <- paste(subdir,'spiec-',as.character(test_number),'.RData',sep="")
save(spiec_record, file=spiec_file)

rm(spiec_glasso, spiec_mb)
rm(spiec_record)
print('after spiec:')
gc()

############################################
# CCLasso
############################################

source("~/LDM-time/kmLDM/program/cclasso.R")
cclasso_result <- cclasso(x=data, counts=T)
cclasso_cov <- cclasso_result$cov.w
cclasso_cor <- cclasso_result$cor.w
print("cclasso result:")
print(cclasso_cov)

cclasso_record <- list(cclasso_cov, cclasso_cor, cclasso_result)

cclasso_file <- paste(subdir,'cclasso-',as.character(test_number),'.RData',sep="")
save(cclasso_record, file=cclasso_file)

rm(cclasso_record, cclasso_result, cclasso_cov, cclasso_cor)
print('after cclasso:')
gc()

# load file
load(pcc_file)
load(scc_file)
load(gl_file)
load(sparcc_file)
load(ccrepe_file)
load(cclasso_file)
load(spiec_file)

otherResults <- list(pcc_record, spearman_record, gl_record, ml_record, sparcc_record, ccrepe_record, spiec_record, cclasso_record)

file.remove(pcc_file, scc_file, gl_file, sparcc_file, ccrepe_file, cclasso_file, spiec_file)

  return(otherResults)
}
