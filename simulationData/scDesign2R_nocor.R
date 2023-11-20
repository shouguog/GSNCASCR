##R4.1.2
rm(list = ls())
library(scDesign2)
library(copula)    # corKendall
library(Rtsne)
library(plyr)      # mapvalues
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggplot2);
theme_set(theme_bw())

# load data -----------------------------------------------------------------------------
data_mat <- readRDS(system.file("extdata", "mouse_sie_10x.rds", package = "scDesign2"))
# remove spike-in -----------------------------------------------------------------------
nonspikes <- which(!grepl("ercc", rownames(data_mat), ignore.case = TRUE))
print(paste("number of spike-ins:", nrow(data_mat)-length(nonspikes)))
#> [1] "number of spike-ins: 9"
data_mat <- data_mat[nonspikes, ,drop = FALSE]
# explore basic structure of data -------------------------------------------------------
dim(data_mat)
#> [1] 15962  7216
table(colnames(data_mat))

unique_cell_type <- names(table(colnames(data_mat)))
set.seed(1)
train_idx <- unlist(sapply(unique_cell_type, function(x){
  cell_type_idx <- which(colnames(data_mat) == x)
  n_cell_total <- length(cell_type_idx)
  sample(cell_type_idx, floor(n_cell_total/2))
}))
traincount <- data_mat[, train_idx]
testcount <- data_mat[, -train_idx]

save(traincount, testcount, file = "ind_train_test.RData")

# set function parameter values ---------------------------------------------------------
n_cell_new <- ncol(testcount[, colnames(testcount) == 'Stem'])

# fit model and simulate data -----------------------------------------------------------
set.seed(1)
ind_result <- fit_model_scDesign2(traincount, 'Stem', sim_method = 'ind')
sim_count_ind <- simulate_count_scDesign2(ind_result, n_cell_new, sim_method = 'ind')
# save the model parameters and the simulated data --------------------------------------
saveRDS(ind_result, file = 'ind_result_Stem_demo.rds')
saveRDS(sim_count_ind, file = 'sim_count_ind_Stem_demo.rds')

##To see how well the simulated data can mimic real data, we will do the following two types of evaluations: 
##(1) Evaluation of the fitting of the marginal distributions, and (2) Evaluation of the characterization of the 
##gene correlations. To evaluate the fitting of the marginal distributions, we will calculate some key gene-wise 
##and cell-wise summary statistics, and compare the distributions of these summary statistics between real data and simulated data.

# a function for computing the marginal stats -------------------------------------------
get_stats <- function(mat, group, log_trans = TRUE){
  mean <- rowMeans(mat)
  var <- apply(mat,1,var)
  cv <- sqrt(var)/mean
  zero_gene <- rowSums(mat < 1e-5)/ncol(mat)
  zero_cell <- colSums(mat < 1e-5)/nrow(mat)
  libsize <- colSums(mat)
  
  if(log_trans){
    mean <- log10(mean + 1)
    var <- log10(var + 1)
    libsize <- log10(libsize + 1)
  }
  summs <- list(mean = mean, var = var, cv = cv, drop_gene = zero_gene,
                drop_cell = zero_cell, libsize = libsize)
  summs = lapply(1:length(summs), function(i){
    data.frame(value = summs[[i]], measure = names(summs)[i], group = group,
               stringsAsFactors = FALSE)
  })
  summs = Reduce(rbind, summs)
  return(summs)
}

# subset traincount and testcount to include only the selected cell type ----------------
traincount_sel <- traincount[, colnames(traincount) == 'Stem']
testcount_sel <- testcount[, colnames(testcount) == 'Stem']
# compute the marginal stats ------------------------------------------------------------
stats_train <- get_stats(traincount_sel, 'training')
stats_test <- get_stats(testcount_sel, 'test')
stats_scDesign2 <- get_stats(sim_count_ind, 'scDesign2')
# organize the stat values as input for ggplot2 -----------------------------------------
stats_dat <- rbind(stats_train, stats_test, stats_scDesign2)
stats_dat$group <- factor(stats_dat$group, levels = c('training', 'test', 'scDesign2'))
measures1 <-  c("mean", "var", "cv", "drop_gene",
                "drop_cell", "libsize")
measures2 <-  c("gene mean", "gene variance", "gene cv",
                "gene zero prop.", "cell zero prop.", "cell library size")
stats_dat$measure <- factor(stats_dat$measure, levels = measures1)
stats_dat$measure <- mapvalues(stats_dat$measure, from = measures1, to = measures2)
# create violin-plots to compare the marginal stat values -------------------------------
stats_plot <- ggplot(stats_dat, aes(x = group, y = value)) +
  geom_violin(scale = 'width', trim = TRUE) +
  facet_wrap(~measure, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("")
print(stats_plot)

##From the above violin-plots, we can see that the distributions of the key summary statistics of the simulated data by scDesign2 
##are very similar to those of the real data (the test data and the training data).

##Next, we will evaluate how well scDesign2 can capture gene correlations. To do this, we will compare the correlation heatmaps 
##among the highly expressed genes between the simulated data and real data. Two types of correlations are considered, Pearson 
##correlation and Kendall’s tau. Pearson correlation is used since it is the most commonly used type of correlation. 
##Kendall’s tau is chosen because it is a rank based correlation that accounts for ties, which is more suitable for count data.

# select the top 100 highly expressed genes ---------------------------------------------
gene_mean <- apply(testcount_sel, 1, mean)
cutoff <- 100
gene_sel <- order(gene_mean, decreasing = TRUE)[1:cutoff]

# two functions for calculating the correlation matrix(-ces) of selected genes ----------
get_cor_mat <- function(x, cor_fun){
  sub_mat <- x[gene_sel, ]
  cor_fun(t(sub_mat))
}
get_heatmap_dat <- function(mat_list, cor_fun){
  cor_mat_list <- lapply(mat_list, get_cor_mat, cor_fun)
  # reorder cor_mat entries according to hierarchical clustering result
  cor_mat_list <- lapply(cor_mat_list, function(x){
    x[hclust_result$order, hclust_result$order]})
  # organize the cor values as input for ggplot2
  cor_melted <- lapply(cor_mat_list, melt)
  cor_dat <- Reduce(rbind, cor_melted)
  cor_dat$group <- unlist(lapply(1:length(group_list), function(x){
    rep(group_list[[x]], nrow(cor_melted[[x]]))
  }))
  return(cor_dat)
}

# calculate the correlations and organize as input for ggplot2 --------------------------
rownames(sim_count_ind) <- rownames(traincount)
mat_list <- list(train = traincount_sel, test = testcount_sel, scDesign2 = sim_count_ind)
hclust_result <- hclust(as.dist(1-get_cor_mat(mat_list$test, cor)))
group_list <- c('training data', 'test data', 'scDesign2')

cor_dat <- get_heatmap_dat(mat_list, cor)
tau_dat <- get_heatmap_dat(mat_list, corKendall)

cor_tau_dat <- rbind(cor_dat, tau_dat)
cor_tau_dat$group <- factor(cor_tau_dat$group, levels = group_list)
cor_tau_dat$cor_type <- factor(c(rep('Pearson\nCorrelation', nrow(cor_dat)),
                                 rep('Kendall\'s\ntau', nrow(tau_dat))),
                               levels = c('Pearson\nCorrelation', 'Kendall\'s\ntau'))

# create heatmaps to display the correlation values -------------------------------------
cor_tau_plot <- ggplot(cor_tau_dat, aes(Var2, Var1, fill = value))+
  facet_grid(vars(cor_type), vars(group)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",name="") +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size=15)) +
  xlab("") + ylab("") + coord_fixed()
png("ind_cor.png", width = 3000, height = 2000, res = 200)
print(cor_tau_plot)
dev.off()


