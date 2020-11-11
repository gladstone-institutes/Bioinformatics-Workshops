##With only the specification of the Coef and not contrast. Clusters have to be specific in terms of character numbers
pbDS_update <- function (pb, method = c("edgeR", "DESeq2", "limma-trend", "limma-voom"), 
                         design = design, coef = Coef, contrast = NULL, min_cells = 10, 
                         verbose = TRUE) 
{
  require(edgeR)
  # method <- match.arg(method)
  # .check_pbs(pb, check_by = TRUE)
  # .check_args_pbDS(as.list(environment()))
  # if (is.null(design)) {
  #   formula <- ~group_id
  #   cd <- as.data.frame(colData(pb))
  #   design <- model.matrix(formula, cd)
  #   colnames(design) <- levels(pb$group_id)
  # }
  # if (is.null(coef) & is.null(contrast)) {
  #   c <- colnames(design)[ncol(design)]
  #   contrast <- makeContrasts(contrasts = c, levels = design)
  # }
  # if (!is.null(contrast)) {
  #   coef <- NULL
  #   ct <- "contrast"
  #   names(cs) <- cs <- colnames(contrast)
  # }
  # else if (!is.null(coef)) {
  #   ct <- "coef"
  #   cs <- vapply(coef, function(i) paste(colnames(design)[i], 
  #                                        collapse = "-"), character(1))
  #   names(cs) <- names(coef) <- cs
  # }
  # print(coef)
  n_cells <- metadata(pb)$n_cells
  kids <- assayNames(pb)
  names(kids) <- kids
  res <- lapply(kids, function(k) {
    if (verbose) 
      cat(k, "..", sep = "")
    rmv <- n_cells[k, ] < min_cells
    y <- assays(pb)[[k]][, !rmv]
    d <- design[!rmv,]
    if (method == "DESeq2") {
      mode(y) <- "integer"
      cd <- colData(pb)[!rmv, , drop = FALSE]
      y <- DESeqDataSetFromMatrix(y, cd, d)
      y <- suppressMessages(DESeq(y))
      tt <- lapply(cs, function(c) {
        res <- results(y, contrast = contrast[, c])
        .res_df(k, res, ct, c) %>% rename(logFC = "log2FoldChange", 
                                          p_val = "pvalue", p_adj.loc = "padj")
      })
    }
    else {
      if (is.null(d) |  any(colSums(d) < 2))
        #if (!is.null(d) & any(colSums(d) < 2))
        return(NULL)
      if (method == "edgeR") {
        if((rowSums(metadata(pb)$n_cells)==0)[(as.integer(k)+1)]){
          return(NULL)
        }
        else {
          y <- suppressMessages(DGEList(y, group = pb$group_id[!rmv], 
                                        remove.zeros = TRUE))
          
          y <- calcNormFactors(y)
          y <- estimateDisp(y, d)
          fit <- glmQLFit(y, d)
          #      tt <- lapply(cs, function(c) {
          qlf <- glmQLFTest(fit, coef)
          tt <- topTags(qlf, n = Inf, sort.by = "none") %>%
            as.data.frame()
          tt <- data.frame(cluster_id=k, tt)
          colnames(tt)[(ncol(tt)-1):ncol(tt)] <- c("p_val","p_adj.loc") 
          #       })
        }
      }
      else {
        if (method == "limma-trend") {
          trend <- robust <- TRUE
        }
        else if (method == "limma-voom") {
          trend <- robust <- FALSE
          y <- suppressMessages(DGEList(y, remove.zeros = TRUE))
          y <- calcNormFactors(y)
          y <- voom(y, d)
        }
        w <- n_cells[k, !rmv]
        fit <- lmFit(y, d, weights = w)
        tt <- lapply(cs, function(c) {
          cfit <- contrasts.fit(fit, contrast[, c], coef[[c]])
          efit <- eBayes(cfit, trend = trend, robust = robust)
          tt <- topTable(efit, number = Inf, sort.by = "none")
          .res_df(k, tt, ct, c) %>% rename(p_val = "P.Value", 
                                           p_adj.loc = "adj.P.Val")
        })
      }
    }
    return(list(tt = tt, data = y))
  })
  skipped <- vapply(res, is.null, logical(1))
  if (any(skipped) & verbose) 
    message("Cluster(s) ", paste(dQuote(kids[skipped]), collapse = ", "), 
            " skipped due to an insufficient number of cells", 
            " in at least 2 samples per group.")
  res <- res[!skipped]
  kids <- kids[names(res)]
  tt <- map(res, "tt")
  # if (!is.null(cs)) {
  #   tt <- lapply(cs, map, .x = tt)
  #   tt <- .p_adj_global(tt)
  # }
  # else {
  #  tt <- lapply(tt, function(u) add_column(u, .after = "p_adj.loc", 
  #                                         p_adj.glb = p.adjust(u$p_adj.loc)))
  # }
  list(table = tt, data = map(res, "data"), method = method, 
       design = design, contrast = contrast, coef = coef)
}

