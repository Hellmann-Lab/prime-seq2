### This is a collection of custom functions
### general functions

## remove Gencode Gene name extension

remove_Geneversion<-function(countmatrix){
rownames(countmatrix) <-str_split_fixed(rownames(countmatrix), pattern="[.]",n=2)[,1]

if (length(unique((rownames(countmatrix))))!=length(rownames(countmatrix))){
  
  countmatrix<-countmatrix %>% 
    as.data.frame() %>% 
    rownames_to_column(var="Gene_ID") %>% 
    mutate(Gene_ID=str_split_fixed(Gene_ID, pattern="[.]",n=2)[,1]) %>% 
    group_by(Gene_ID) %>% 
    dplyr::summarize(across(where(is.double),sum)) %>% 
    ungroup() %>% 
    column_to_rownames(var="Gene_ID") %>% 
    as.matrix()
  return(countmatrix) 
}else {
 return(countmatrix) 
}

}



#### PowsimR helper functions to extract results for plotting

getEvalDE_df_marginal<-function(evalres,method=c("prime-seq","tru-seq"),count_type=c("exon","inex")){
  dat.marginal <- evalres[grep('*R.marginal', names(evalres))]
  
  names(dat.marginal) <- substr(x = names(dat.marginal), start = 1, stop = 3)
  dat.marginal <- lapply(dat.marginal, "rownames<-", paste0(evalres[['n1']], " vs ", evalres[['n2']]))
  dat.marginal.long <- reshape2::melt(dat.marginal)
  dat.marginal.long$L1 <- factor(dat.marginal.long$L1,
                                 levels = c("TPR", "FNR", "FPR", "TNR", "FDR"))
  
  dat.marginal.long$method<-method
  
  dat.marginal.long$count_type<-count_type
  return(dat.marginal.long)
}


## conditional rates
getEvalDE_list_conditional<-function(evalRes,method=c("prime-seq","tru-seq"),count_type=c("exon","inex")){
  plot_list<-list(dat.genes.calc=NULL,dat.stratified.long=NULL,stratum.name=NULL)
  plot_list$stratum.name <- dplyr::case_when(evalRes$stratify.by == "mean" ~ "(Log2 Mean Expression)",
                                   evalRes$stratify.by == "dispersion" ~ "(Log2 Dispersion)",
                                   evalRes$stratify.by == "lfc" ~ "(Log2 Fold Change)",
                                   evalRes$stratify.by == "lfc_abs" ~ "absolute (Log2 Fold Change)",
                                   evalRes$stratify.by == "dropout" ~ "(Gene Dropout Rate)")
  # strata genes
  strata <- evalRes$strata.levels
  N <- length(evalRes$n1)
  dat.genes <- list("Ngenes"=evalRes$stratagenes[,N,],
                    'DEgenes'=evalRes$stratadiffgenes[,N,])
  dat.genes <- lapply(dat.genes, "rownames<-", strata)
  dat.genes.long <- reshape2::melt(dat.genes)
  dat.genes.calc <- dat.genes.long %>%
    dplyr::group_by(.data$Var1, .data$L1) %>%
    dplyr::summarise(Expectation=mean(.data$value),
                     Deviation=stats::sd(.data$value),
                     Error=stats::sd(.data$value)/sqrt(dplyr::n())) %>%
    dplyr::ungroup()
  dat.genes.calc$method<-method
  dat.genes.calc$count_type<-count_type
  plot_list$dat.genes.calc<-dat.genes.calc
  refval <- data.frame(L1 = c("FDR", "TPR"),
                       ref = c(evalRes$alpha.nominal, 0.8))
  
  
  dat.stratified <- evalRes[grep('*R$', names(evalRes))]
  
  dat.stratified <- lapply(dat.stratified, "dimnames<-",
                           list(strata, paste0(evalRes[['n1']], " vs ", evalRes[['n2']]),
                                NULL))
  dat.stratified.long <- reshape2::melt(dat.stratified)
  dat.stratified.long$L1 <- factor(dat.stratified.long$L1,
                                   levels = c("TPR", "FNR", "FPR", "TNR", "FDR"))
  
  dat.stratified.long$method<-method
  
  dat.stratified.long$count_type<-count_type
  
  plot_list$dat.stratified.long<-dat.stratified.long
  return(plot_list)
}

### Get Mean dispersion plot from PowsimR estparameters

Get_Mean_Disp<- function(estParamRes){
  ..density.. = NULL
  
  plot_list<-list(gene_id=NULL,meanvsdisp.dat=NULL,meanvsdisp.fdat=NULL,cdisp=NULL)
  plot_list$gene_id<-estParamRes[["Fit"]][["Filtered"]][["sharedgenes"]]
  plot_list$meanvsdisp.dat <- data.frame(Mean=estParamRes$Fit$Filtered$meandispfit$model$x[,"x"],
                                         Dispersion=estParamRes$Fit$Filtered$meandispfit$model$y)
  plot_list$meanvsdisp.fdat <- data.frame(Mean=estParamRes$Fit$Filtered$meandispfit$x,
                                          Dispersion=estParamRes$Fit$Filtered$meandispfit$y,
                                          Upper=estParamRes$Fit$Filtered$meandispfit$upper,
                                          Lower=estParamRes$Fit$Filtered$meandispfit$lower)
  
  plot_list$cdisp <- log2(estParamRes$Parameters$Filtered$common.dispersion+1)
  
  return(plot_list)
}







evaluateDE2<-function (simRes, alpha.type = c("adjusted", "raw"), MTC = c("BY", 
                                                                          "BH", "holm", "hochberg", "hommel", "bonferroni", "Storey", 
                                                                          "IHW"), alpha.nominal = 0.1, stratify.by = c("mean", "dispersion", 
                                                                                                                       "dropout", "lfc","lfc_abs"), strata.probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 
                                                                                                                                                                     0.6, 0.7, 0.8, 0.9), filter.by = c("none", "mean", "dispersion", 
                                                                                                                                                                                                        "dropout"), strata.filtered = 1, target.by = c("lfc", "effectsize"), 
                       delta = 0, Table = TRUE) 
{
  alpha.type = match.arg(alpha.type)
  MTC = match.arg(MTC)
  stratify.by = match.arg(stratify.by)
  filter.by = match.arg(filter.by)
  target.by = match.arg(target.by)
  Nreps1 = simRes$SimSetup$n1
  Nreps2 = simRes$SimSetup$n2
  ngenes = simRes$DESetup$ngenes
  DEids = simRes$DESetup$DEid
  lfcs = simRes$DESetup$pLFC
  tlfcs = lapply(1:length(lfcs), function(i) {
    lfcs[[i]]
  })
  nsims = simRes$DESetup$nsims
  estmeans = simRes$estParamRes$Parameters[[simRes$DESetup$Draw$MoM]]$means
  estdisps = simRes$estParamRes$Parameters[[simRes$DESetup$Draw$MoM]]$dispersion
  estdropout = simRes$estParamRes$Parameters[[simRes$DESetup$Draw$MoM]]$gene.dropout
  mu = simRes$SimulateRes$mu
  disp = simRes$SimulateRes$disp
  dropout = simRes$SimulateRes$dropout
  elfc = simRes$SimulateRes$elfc
  DEmethod = simRes$Pipeline$DEmethod
  pvalue = simRes$SimulateRes$pvalue
  fdr = simRes$SimulateRes$fdr
  tmp.ecdf.mean = stats::ecdf(log2(estmeans + 1))
  tmp.quantile.mean = stats::quantile(tmp.ecdf.mean, probs = strata.probs)
  strata.mean = unique(c(0, unname(tmp.quantile.mean), Inf))
  strata.mean = unique(round(strata.mean, digits = 2))
  tmp.ecdf.disps = stats::ecdf(log2(estdisps + 1))
  tmp.quantile.disps = stats::quantile(tmp.ecdf.disps, probs = strata.probs)
  strata.disps = unique(c(0, unname(tmp.quantile.disps), Inf))
  strata.disps = unique(round(strata.disps, digits = 2))
  tmp.ecdf.drop = stats::ecdf(estdropout)
  tmp.quantile.drop = stats::quantile(tmp.ecdf.drop, probs = strata.probs)
  strata.drop = unique(c(0, unname(tmp.quantile.drop), 1))
  strata.drop = unique(round(strata.drop, digits = 2))
  tmp.ecdf.lfc = stats::ecdf(unique(unlist(tlfcs)))
  tmp.quantile.lfc = stats::quantile(tmp.ecdf.lfc, probs = strata.probs)
  strata.lfc = unique(c(-Inf, unname(tmp.quantile.lfc), Inf))
  strata.lfc = unique(round(strata.lfc, digits = 2))
  tmp.ecdf.lfc.abs = stats::ecdf(unique(abs(unlist(tlfcs))))
  tmp.quantile.lfc.abs = stats::quantile(tmp.ecdf.lfc.abs, probs = strata.probs)
  strata.lfc.abs = unique(c(0, unname(tmp.quantile.lfc.abs), Inf))
  strata.lfc.abs = unique(round(strata.lfc.abs, digits = 2))
  
  if (stratify.by == "mean") {
    nr = length(strata.mean) - 1
  }
  if (stratify.by == "dispersion") {
    nr = length(strata.disps) - 1
  }
  if (stratify.by == "dropout") {
    nr = length(strata.drop) - 1
  }
  if (stratify.by == "lfc") {
    nr = length(strata.lfc) - 1
  }
  if (stratify.by == "lfc_abs") {
    nr = length(strata.lfc.abs) - 1
  }
  if (filter.by %in% c("mean", "dispersion", "dropout")) {
    if (filter.by == stratify.by) {
      nstrata = nr - strata.filtered
    }
    if (!filter.by == stratify.by) {
      nstrata = nr
    }
  }
  if (filter.by == "none") {
    nstrata = nr
  }
  TP = TN = FP = FN = TPR = TNR = FPR = FNR = FDR = xgrl = xgrld = array(NA, 
                                                                         dim = c(nstrata, length(Nreps1), nsims))
  TP.marginal = TN.marginal = FP.marginal = FN.marginal = TPR.marginal = TNR.marginal = FPR.marginal = FNR.marginal = FDR.marginal = matrix(NA, 
                                                                                                                                            length(Nreps1), nsims)
  for (i in 1:nsims) {
    for (j in seq(along = Nreps1)) {
      Nrep1 = Nreps1[j]
      Nrep2 = Nreps2[j]
      DEid = DEids[[i]]
      lfc = lfcs[[i]]
      lfc_abs = abs(lfcs[[i]])
      Zg = Zg2 = rep(0, ngenes)
      Zg[DEid] = 1
      if (delta == 0) {
        Zg2 = Zg
      }
      if (!delta == 0) {
        if (target.by == "lfc") {
          ix = abs(lfc) > delta
        }
        else if (target.by == "effectsize") {
          effectsize = lfc/sqrt(1/((log2(mu[, j, i] + 
                                           1) + log2(disp[, j, i] + 1))))
          ix = abs(effectsize) > delta
        }
        Zg2[ix] = 1
      }
      X.bar1 = mu[, j, i]
      ix.keep.mean = which(!is.na(X.bar1))
      xgr.mean = cut(log2(X.bar1[ix.keep.mean] + 1), strata.mean)
      xgrd.mean = cut(log2(X.bar1[DEid] + 1), strata.mean)
      X.disp1 = disp[, j, i]
      ix.keep.disps = which(!is.na(X.disp1))
      xgr.disps = cut(log2(X.disp1[ix.keep.disps] + 1), 
                      strata.disps)
      xgrd.disps = cut(log2(X.disp1[DEid] + 1), strata.disps)
      X.drop1 = dropout[, j, i]
      ix.keep.drop = which(!is.na(X.drop1))
      xgr.drop = cut(X.drop1[ix.keep.drop], strata.drop)
      xgrd.drop = cut(X.drop1[DEid], strata.drop)
      X.lfc1 = elfc[, j, i]
      ix.keep.lfc = which(!is.na(X.lfc1))
      xgr.lfc = cut(X.lfc1[ix.keep.lfc], strata.lfc)
      xgrd.lfc = cut(X.lfc1[DEid], strata.lfc)
      X.lfc_abs1 = abs(elfc[, j, i])
      ix.keep.lfc_abs = which(!is.na(X.lfc_abs1))
      xgr.lfc_abs = cut(X.lfc_abs1[ix.keep.lfc_abs], strata.lfc.abs)
      xgrd.lfc_abs = cut(X.lfc_abs1[DEid], strata.lfc.abs)
      if (stratify.by == "mean") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.mean = ix.keep.mean[!(xgr.mean %in% 
                                           lev.mean[strata.filt.mean])]
          fix.dekeep.mean = intersect(fix.keep.mean, 
                                      DEid)
          fxgr.mean = cut(log2(X.bar1[fix.keep.mean] + 
                                 1), strata.mean[-strata.filt.mean])
          fxgrd.mean = cut(log2(X.bar1[fix.dekeep.mean] + 
                                  1), strata.mean[-strata.filt.mean])
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.mean = ix.keep.mean[!(xgr.disps %in% 
                                           lev.disps[strata.filt.disps])]
          fix.dekeep.mean = intersect(fix.keep.mean, 
                                      DEid)
          fxgr.mean = cut(log2(X.bar1[fix.keep.mean] + 
                                 1), strata.mean)
          fxgrd.mean = cut(log2(X.bar1[fix.dekeep.mean] + 
                                  1), strata.mean)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.mean = ix.keep.mean[!(xgr.drop %in% 
                                           lev.drop[strata.filt.drop])]
          fix.dekeep.mean = intersect(fix.keep.mean, 
                                      DEid)
          fxgr.mean = cut(log2(X.bar1[fix.keep.mean] + 
                                 1), strata.mean)
          fxgrd.mean = cut(log2(X.bar1[fix.dekeep.mean] + 
                                  1), strata.mean)
        }
        if (filter.by == "none") {
          fix.keep.mean = ix.keep.mean
          fxgr.mean = xgr.mean
          fxgrd.mean = xgrd.mean
        }
      }
      if (stratify.by == "dispersion") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.disps = ix.keep.disps[!(xgr.mean %in% 
                                             lev.mean[strata.filt.mean])]
          fix.dekeep.disps = intersect(fix.keep.disps, 
                                       DEid)
          fxgr.disps = cut(log2(X.disp1[fix.keep.disps] + 
                                  1), strata.disps)
          fxgrd.disps = cut(log2(X.disp1[fix.dekeep.disps] + 
                                   1), strata.disps)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          lev.filt.disps = c((length(lev.disps) - strata.filtered + 
                                1):length(lev.disps))
          fix.keep.disps = ix.keep.disps[!(xgr.disps %in% 
                                             lev.disps[lev.filt.disps])]
          fix.dekeep.disps = intersect(fix.keep.disps, 
                                       DEid)
          strata.filt.disps = length(strata.disps) - 
            strata.filtered
          fxgr.disps = cut(log2(X.disp1[fix.keep.disps] + 
                                  1), strata.disps[1:strata.filt.disps])
          fxgrd.disps = cut(log2(X.disp1[fix.dekeep.disps] + 
                                   1), strata.disps[1:strata.filt.disps])
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.disps = ix.keep.disps[!(xgr.drop %in% 
                                             lev.drop[strata.filt.drop])]
          fix.dekeep.disps = intersect(fix.keep.disps, 
                                       DEid)
          fxgr.disps = cut(X.disp1[fix.keep.disps], strata.disps)
          fxgrd.disps = cut(X.disp1[fix.dekeep.disps], 
                            strata.disps)
        }
        if (filter.by == "none") {
          fix.keep.disps = ix.keep.disps
          fxgr.disps = xgr.disps
          fxgrd.disps = xgrd.disps
        }
      }
      if (stratify.by == "dropout") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.drop = ix.keep.drop[!(xgr.mean %in% 
                                           lev.mean[strata.filt.mean])]
          fix.dekeep.drop = intersect(fix.keep.drop, 
                                      DEid)
          fxgr.drop = cut(X.drop1[fix.keep.drop], strata.drop)
          fxgrd.drop = cut(X.drop1[fix.dekeep.drop], 
                           strata.drop)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.drop = ix.keep.drop[!(xgr.disps %in% 
                                           lev.disps[strata.filt.disps])]
          fix.dekeep.drop = intersect(fix.keep.drop, 
                                      DEid)
          fxgr.drop = cut(X.drop1[fix.keep.drop], strata.drop)
          fxgrd.drop = cut(X.drop1[fix.dekeep.drop], 
                           strata.drop)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.drop = ix.keep.drop[!(xgr.drop %in% 
                                           lev.drop[strata.filt.drop])]
          fix.dekeep.drop = intersect(fix.keep.drop, 
                                      DEid)
          strata.filt.drop = length(strata.drop) - strata.filtered
          fxgr.drop = cut(X.drop1[fix.keep.drop], strata.drop[1:strata.filt.drop])
          fxgrd.drop = cut(X.drop1[fix.dekeep.drop], 
                           strata.drop[1:strata.filt.drop])
        }
        if (filter.by == "none") {
          fix.keep.drop = ix.keep.drop
          fxgr.drop = xgr.drop
          fxgrd.drop = xgrd.drop
        }
      }
      if (stratify.by == "lfc") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.lfc = ix.keep.lfc[!(xgr.mean %in% 
                                         lev.mean[strata.filt.mean])]
          fix.dekeep.lfc = intersect(fix.keep.lfc, DEid)
          fxgr.lfc = cut(X.lfc1[fix.keep.lfc], strata.lfc)
          fxgrd.lfc = cut(X.lfc1[fix.dekeep.lfc], strata.lfc)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.lfc = ix.keep.lfc[!(xgr.disps %in% 
                                         lev.mean[strata.filt.disps])]
          fix.dekeep.lfc = intersect(fix.keep.lfc, DEid)
          fxgr.lfc = cut(X.lfc1[fix.keep.lfc], strata.lfc)
          fxgrd.lfc = cut(X.lfc1[fix.dekeep.lfc], strata.lfc)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.lfc = ix.keep.lfc[!(xgr.drop %in% 
                                         lev.drop[strata.filt.drop])]
          fix.dekeep.lfc = intersect(fix.keep.lfc, DEid)
          fxgr.lfc = cut(X.lfc1[fix.keep.lfc], strata.lfc)
          fxgrd.lfc = cut(X.lfc1[fix.dekeep.lfc], strata.lfc)
        }
        if (filter.by == "none") {
          fix.keep.lfc = ix.keep.lfc
          fxgr.lfc = xgr.lfc
          fxgr.lfc = xgr.lfc
          fxgrd.lfc = xgrd.lfc
        }
      }
      if (stratify.by == "lfc_abs") {
        if (filter.by == "mean") {
          lev.mean = levels(xgr.mean)
          strata.filt.mean = c(1:strata.filtered)
          fix.keep.lfc_abs = ix.keep.lfc_abs[!(xgr.mean %in% 
                                                 lev.mean[strata.filt.mean])]
          fix.dekeep.lfc_abs = intersect(fix.keep.lfc_abs, DEid)
          fxgr.lfc_abs = cut(X.lfc_abs1[fix.keep.lfc_abs], strata.lfc.abs)
          fxgrd.lfc_abs = cut(X.lfc_abs1[fix.dekeep.lfc_abs], strata.lfc.abs)
        }
        if (filter.by == "dispersion") {
          lev.disps = levels(xgr.disps)
          strata.filt.disps = c((length(lev.disps) - 
                                   strata.filtered + 1):length(lev.disps))
          fix.keep.lfc_abs = ix.keep.lfc_abs[!(xgr.disps %in% 
                                                 lev.mean[strata.filt.disps])]
          fix.dekeep.lfc_abs = intersect(fix.keep.lfc_abs, DEid)
          fxgr.lfc_abs = cut(X.lfc_abs1[fix.keep.lfc_abs], strata.lfc.abs)
          fxgrd.lfc_abs = cut(X.lfc_abs1[fix.dekeep.lfc_abs], strata.lfc.abs)
        }
        if (filter.by == "dropout") {
          lev.drop = levels(xgr.drop)
          strata.filt.drop = c((length(lev.drop) - strata.filtered + 
                                  1):length(lev.drop))
          fix.keep.lfc_abs = ix.keep.lfc_abs[!(xgr.drop %in% 
                                                 lev.drop[strata.filt.drop])]
          fix.dekeep.lfc_abs = intersect(fix.keep.lfc_abs, DEid)
          fxgr.lfc_abs = cut(X.lfc_abs1[fix.keep.lfc_abs], strata.lfc.abs)
          fxgrd.lfc_abs = cut(X.lfc_abs1[fix.dekeep.lfc_abs], strata.lfc.abs)
        }
        if (filter.by == "none") {
          fix.keep.lfc_abs = ix.keep.lfc_abs
          fxgr.lfc_abs = xgr.lfc_abs
          fxgr.lfc_abs = xgr.lfc_abs
          fxgrd.lfc_abs = xgrd.lfc_abs
        }
      }
      if (stratify.by == "mean") {
        strata = strata.mean
        xgr = fxgr.mean
        xgrd = fxgrd.mean
        ix.keep = fix.keep.mean
      }
      if (stratify.by == "dispersion") {
        strata = strata.disps
        xgr = fxgr.disps
        xgrd = fxgrd.disps
        ix.keep = fix.keep.disps
      }
      if (stratify.by == "dropout") {
        strata = strata.drop
        xgr = fxgr.drop
        xgrd = fxgrd.drop
        ix.keep = fix.keep.drop
      }
      if (stratify.by == "lfc") {
        strata = strata.lfc
        xgr = fxgr.lfc
        xgrd = fxgrd.lfc
        ix.keep = fix.keep.lfc
      }
      if (stratify.by == "lfc_abs") {
        strata = strata.lfc.abs
        xgr = fxgr.lfc_abs
        xgrd = fxgrd.lfc_abs
        ix.keep = fix.keep.lfc_abs
      }
      if (alpha.type == "raw") {
        if (DEmethod %in% c("edgeR-QL", "edgeR-LRT", 
                            "limma-voom", "limma-trend", "NBPSeq", "T-Test", 
                            "DESeq2", "ROTS", "MAST", "scde", "BPSC", "scDD", 
                            "monocle", "DECENT", "edgeR-zingeR", "edgeR-ZINB-WaVE", 
                            "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          x = pvalue[ix.keep, j, i]
          x[is.na(x)] = 1
        }
        if (DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod, 
                         " only provides adjusted p-values."))
          x = fdr[ix.keep, j, i]
          x[is.na(x)] = 1
        }
      }
      if (alpha.type == "adjusted") {
        if (DEmethod %in% c("edgeR-QL", "edgeR-LRT", 
                            "limma-voom", "limma-trend", "NBPSeq", "T-Test", 
                            "DESeq2", "ROTS", "MAST", "scde", "BPSC", "scDD", 
                            "monocle", "DECENT", "edgeR-zingeR", "edgeR-ZINB-WaVE", 
                            "DESeq2-zingeR", "DESeq2-ZINB-WaVE")) {
          pval = pvalue[ix.keep, j, i]
          meanexpr = log2(mu[ix.keep, j, i] + 1)
          if (MTC %in% stats::p.adjust.methods) {
            x = stats::p.adjust(pval, method = MTC)
            x[is.na(x)] = 1
          }
          if (MTC %in% "Storey") {
            tmp.p = pval[!is.na(pval)]
            tmp.q = qvalue::qvalue(p = tmp.p)$qvalues
            x = rep(NA, length(pval))
            x[!is.na(pval)] = tmp.q
            x[is.na(x)] = 1
          }
          if (MTC %in% "IHW") {
            in.dat = data.frame(pvalue = pval, meanexpr = meanexpr)
            tmp = IHW::ihw(pvalue ~ meanexpr, data = in.dat, 
                           alpha = alpha.nominal)
            x = IHW::adj_pvalues(tmp)
            x[is.na(x)] = 1
          }
        }
        if (DEmethod %in% c("baySeq", "NOISeq", "EBSeq")) {
          message(paste0("The DE method ", DEmethod, 
                         " only provides adjusted p-values."))
          x = fdr[ix.keep, j, i]
          x[is.na(x)] = 1
        }
      }
      Zg = Zg[ix.keep]
      Zg2 = Zg2[ix.keep]
      xgrl[, j, i] = table(xgr)
      xgrld[, j, i] = table(xgrd)
      error.mat = powsimR:::.error.matrix(p = x, p.crit = alpha.nominal, 
                                          Zg = Zg, Zg2 = Zg2, xgr = xgr)
      TP[, j, i] = error.mat$TP
      TN[, j, i] = error.mat$TN
      FP[, j, i] = error.mat$FP
      FN[, j, i] = error.mat$FN
      TP.marginal[j, i] = error.mat$TP.marginal
      TN.marginal[j, i] = error.mat$TN.marginal
      FP.marginal[j, i] = error.mat$FP.marginal
      FN.marginal[j, i] = error.mat$FN.marginal
      TPR[, j, i] = error.mat$TPR
      TNR[, j, i] = error.mat$TNR
      FPR[, j, i] = error.mat$FPR
      FNR[, j, i] = error.mat$FNR
      FDR[, j, i] = error.mat$FDR
      TPR.marginal[j, i] = error.mat$TPR.marginal
      TNR.marginal[j, i] = error.mat$TNR.marginal
      FPR.marginal[j, i] = error.mat$FPR.marginal
      FNR.marginal[j, i] = error.mat$FNR.marginal
      FDR.marginal[j, i] = error.mat$FDR.marginal
    }
  }
  output <- list(stratagenes = xgrl, stratadiffgenes = xgrld, 
                 TN = TN, TP = TP, FP = FP, FN = FN, TN.marginal = TN.marginal, 
                 TP.marginal = TP.marginal, FP.marginal = FP.marginal, 
                 FN.marginal = FN.marginal, TNR = TNR, TPR = TPR, FPR = FPR, 
                 FNR = FNR, FDR = FDR, TNR.marginal = TNR.marginal, TPR.marginal = TPR.marginal, 
                 FPR.marginal = FPR.marginal, FNR.marginal = FNR.marginal, 
                 FDR.marginal = FDR.marginal, alpha.type = alpha.type, 
                 MTC = ifelse(alpha.type == "adjusted", MTC, "not applicable"), 
                 alpha.nominal = alpha.nominal, stratify.by = stratify.by, 
                 strata = strata, strata.levels = levels(xgr), target.by = target.by, 
                 n1 = Nreps1, n2 = Nreps2, delta = delta)
  if (Table) {
    printEvalDE(evalRes = output)
  }
  return(output)
}
