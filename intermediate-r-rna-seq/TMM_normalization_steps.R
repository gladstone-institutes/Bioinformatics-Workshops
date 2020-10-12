#The code for the main steps in the following ...
#... are copied from the source code for edgeR::calcNormFactors.

cnts <- y$counts %>% as.matrix()
lib.sizes <- apply(cnts, 2, sum)

cnts_adjst_libsize <- map_dfc(colnames(cnts), 
                              function(x) cnts[, x]/lib.sizes[x]) %>%
  set_colnames(., colnames(cnts))

boxplot(cnts_adjst_libsize)

f <- apply(cnts_adjst_libsize, 2, 
           function(x) quantile(x, p=0.75))

boxplot(cnts_adjst_libsize, ylim = c(0, 10^(-4)))
abline(h = mean(f), lty = "dotted")

ref_sample <- (f - mean(f)) %>% 
  abs() %>% 
  which.min()

#Illustrating normalization of one sample.
illustrate = TRUE
if (illustrate) {
x <- colnames(cnts)[12]

nO <- lib.sizes[x]
nR <- lib.sizes[ref_sample]

obs <- cnts[, x] %>% as.numeric()
ref <- cnts[, ref_sample] %>% as.numeric()

#The M values:
logR <- log2(obs/nO) - log2(ref/nR)         # log ratio of expression, accounting for library size

#The A values:
absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression

#	remove infinite values, cutoff based on A
fin <- is.finite(logR) & is.finite(absE) & (absE > -10^(10))

logR <- logR[fin]
absE <- absE[fin]

print(
ggplot(data.frame(A = absE, 
                  M = logR), 
       aes(A, M)) + 
  geom_point() + 
  geom_smooth() + 
  coord_cartesian(xlim = c(-25, -5), ylim = c(-10.5, 10.5))
)

logratioTrim=.3
sumTrim=0.05

#Remove the genes with the 5% most extreme A values.
#Remove the genes with the 30% most extreme M values
n <- length(logR)
loL <- floor(n * logratioTrim) + 1
hiL <- n + 1 - loL
loS <- floor(n * sumTrim) + 1
hiS <- n + 1 - loS

keep <- (rank(logR)>=loL & rank(logR)<=hiL) & 
  (rank(absE)>=loS & rank(absE)<=hiS)

print(
  ggplot(data.frame(A = absE[keep], 
                  M = logR[keep]), 
       aes(A, M)) + 
  geom_point() + 
  coord_cartesian(xlim = c(-25, -5), ylim = c(-10.5, 10.5)) +
  geom_hline(color = "red", linetype = "dotted", 
           yintercept = mean(logR[keep], na.rm = T))
)

#Normalization factor for the current sample wrt to the reference sample
f <- mean(logR[keep], na.rm=TRUE) 
print("Normalization factor (log2 scale)")
print(f)

f <- 2^f
print("Normalization factor (original scale)")
print(f)

#Compare with the normalization factor calculated by TMM.
print("Normalization factor (from edgeR::calcNormFactors)")
y$samples$norm.factors[12]
}

TMM_norm_factors <- 
  map_dfc(colnames(cnts), 
          function(x, logratioTrim=.3, sumTrim=0.05, 
                   doWeighting=TRUE, Acutoff=-1e10) {
            
            #The following steps are excerpted from edgeR's source.
            nO <- lib.sizes[x]
            nR <- lib.sizes[ref_sample]
            
            obs <- cnts[, x] %>% as.numeric()
            ref <- cnts[, ref_sample] %>% as.numeric()
            
            logR <- log2(obs/nO) - log2(ref/nR)         # log ratio of expression, accounting for library size
            absE <- (log2(obs/nO) + log2(ref/nR))/2  # absolute expression
            v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref   # estimated asymptotic variance
            
            #	remove infinite values, cutoff based on A
            fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
            
            logR <- logR[fin]
            absE <- absE[fin]
            v <- v[fin]
            
            if(max(abs(logR)) < 1e-6) return(1)
            
            #	taken from the original mean() function
            n <- length(logR)
            loL <- floor(n * logratioTrim) + 1
            hiL <- n + 1 - loL
            loS <- floor(n * sumTrim) + 1
            hiS <- n + 1 - loS
            
            keep <- (rank(logR)>=loL & rank(logR)<=hiL) & 
              (rank(absE)>=loS & rank(absE)<=hiS)
            
            if(doWeighting)
              f <- sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE)
            else
              f <- mean(logR[keep], na.rm=TRUE)
            
            #	Results will be missing if the two libraries share no features with positive counts
            #	In this case, return unity
            if(is.na(f)) f <- 0
            2^f
          }) %>%
  data.frame() %>% 
  as.numeric() 

#Rescale norm factors for convenience of interpretation.
rescale <- TMM_norm_factors %>% 
  log() %>%
  mean() %>%
  exp()

TMM_norm_factors %<>% divide_by(., rescale)
TMM_norm_factors

#Compare with the output from calcNormFactors
y$samples$norm.factors

