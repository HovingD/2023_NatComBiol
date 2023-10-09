# Script for the analysis of Snp datasets (elbow and Ig)
# model: value ~ type + timepoint

#renv::activate() # if you want to use renv

# PACKAGES ---------------------------------------------------------------------

library(plyr)             # install.packages('plyr')
library(dplyr)            # install.packages('dplyr')
library(tidyr)            # install.packages('tidyr')
library(magrittr)         # install.packages('magrittr')
library(stringr)          # install.packages('stringr')
library(ggplot2)          # install.packages('ggplot2')
library(writexl)          # install.packages('writexl')
library(lmboot)           # install.packages('lmboot')

# PARAMETERS -------------------------------------------------------------------
                       # Make changes only in this section #
       # in the dataset use lowercase for the file and donor column names #

# path to the folder where the results will be saved
path <- ".//"

# project name
file_name <- '\'

# data set for the analysis
data <- read.csv("./.csv", sep=",")

# regular expression (possible values are: 'elbow', 'Ig')
regex_pattern <- 'elbow'

# column for grouping (ex: 'timepoint', 'group')
group_col <- 'timepoint'

# group_col caregory to be used as baseline for comparisons
group_col_baseline <- 'pre'

# number of bootstrap resamples
B <- 9999

# kind of hypothis testing (possible values are: 'two-tailed', 'upper-tailed', 'lower-tailed')
test_kind <- 'upper-tailed'

# significance level for the test
alpha_level <- 0.05

# kind of p-value (possible values are: 'unadjusted', 'adjusted')
p_kind <- 'adjusted'

# vector including all the cluster names in a desired order
ordered_cluster <- NULL

# vector including all the PS types in a desired order
ordered_type <- NULL

# vector including all the timepoints in a desired order
ordered_timepoint <- c('pre', 'M4')

# TRANSFORM DATA ---------------------------------------------------------------

data %<>% replace_na(list(value=0))

data_long <- data %>%
  dplyr::select(matches(paste0('(^PS.*', regex_pattern, ')|(^file$)|(^donor$)|(^', group_col, '$)'))) %>%
  gather('type', 'value', -c('file', 'donor', eval(group_col))) %>%
  mutate(cluster=str_extract(type, paste0(regex_pattern, '..')),
         type=gsub(paste0('\\.', substr(regex_pattern, 1, 1)), '', str_extract(type, paste0('PS.*\\.', substr(regex_pattern, 1, 1)))),
         file=gsub('^Batch[[:digit:]]+_|\\.csv$', '', file)) %>%
  rename(timepoint=eval(group_col)) %>%
  ddply(~ file + donor + timepoint + type + cluster, summarise, value=mean(value)) %>%
  dplyr::select(-file)

data_ps_filtered <- data_long %>%
  ddply(~ type, summarise, total=sum(value)) %>%
  filter(total > 4) %>% #change to number of minimum cells per PS
  dplyr::select(type) %>%
  unlist() %>% unique()

data_perc <- data_long %>%
  filter(type %in% data_ps_filtered) %>%
  ddply(~ donor + timepoint + type, mutate, value=value / sum(value)) %>%
  replace_na(list(value=0))

type_names <- unique(data_perc$type)
timepoint_names <- unique(data_perc$timepoint)
cluster_names <- unique(data_perc$cluster)

data_changes_outline <-data.frame(kept=paste(unique(data_perc$type), collapse=', '),
                                  excluded=paste(setdiff(unique(data_long$type), unique(data_perc$type)), collapse=', '))

# WILD BOOTSTRAP ---------------------------------------------------------------

set.seed(1234)

k <- ifelse(test_kind == 'two-tailed', 2, 1)

list_wild_boot_type <- lapply(cluster_names, function(j) {
  cluster_data_perc <- subset(data_perc, cluster == j)
  
  fit <- lm(value ~ type + timepoint, data=cluster_data_perc)
  
  n_type <- length(type_names)
  ps_no_intercept <- na.omit(str_extract(names(coef(fit)), 'PS.*'))
  ps_intercept <- setdiff(type_names, ps_no_intercept)
  
  wild_boot <- lapply(c(ps_no_intercept, 'Intercept'), function(i) {
    cat('running wild bootstrap for', j, ifelse(i == 'Intercept', ps_intercept, i), '...\n')

    regex_string <- ifelse(i == 'Intercept', '\\(.*\\)', paste0('type', i, '$|\\(.*\\)'))
    type_ps <- na.omit(str_extract(names(coef(fit)), regex_string))
    type_coefs <- na.omit(c(str_extract(names(coef(fit)), 'type.*|\\(.*\\)'), type_ps))
    type_weights <- c(-c(n_type, rep(1, n_type - 1)) / n_type, rep(1, length(type_ps)))
  
    suppressWarnings({
      fit_wild_boot <- wild.boot(value ~ type + timepoint, data=cluster_data_perc, B=B)
    })
    
    obs_est <- sum(coef(fit)[type_coefs] * type_weights)
    boot_est <- c(fit_wild_boot$bootEstParam[, type_coefs] %*% matrix(type_weights, ncol=1))

    limits <- c(alpha_level / k, 1 - alpha_level / k)
    boot_ci <- quantile(boot_est, probs=limits) %>% matrix(nrow=1)
    boot_ci[1] <- ifelse(test_kind == 'lower-tailed', -Inf, boot_ci[1])
    boot_ci[2] <- ifelse(test_kind == 'upper-tailed', Inf, boot_ci[2])
    colnames(boot_ci) %<>% paste0('ci_limit_', limits)
    colnames(boot_ci)[1] <- ifelse(test_kind == 'lower-tailed', 'ci_limit_minusInf', colnames(boot_ci)[1])
    colnames(boot_ci)[2] <- ifelse(test_kind == 'upper-tailed', 'ci_limit_Inf', colnames(boot_ci)[2])
    signif_ci <- ifelse(prod(sign(boot_ci)) > 0, 1, 0)

    adj <- length(ps_no_intercept) + 1
    limits_bonf <- c((alpha_level / adj) / k, 1 - (alpha_level / adj) / k)
    boot_ci_bonf <- quantile(boot_est, probs=limits_bonf) %>% matrix(nrow=1)
    boot_ci_bonf[1] <- ifelse(test_kind == 'lower-tailed', -Inf, boot_ci_bonf[1])
    boot_ci_bonf[2] <- ifelse(test_kind == 'upper-tailed', Inf, boot_ci_bonf[2])
    colnames(boot_ci_bonf) %<>% paste0('ci.adj_limit_', limits_bonf)
    colnames(boot_ci_bonf)[1] <- ifelse(test_kind == 'lower-tailed', 'ci.adj_limit_minusInf', colnames(boot_ci_bonf)[1])
    colnames(boot_ci_bonf)[2] <- ifelse(test_kind == 'upper-tailed', 'ci.adj_limit_Inf', colnames(boot_ci_bonf)[2])
    signif_ci_adj <- ifelse(prod(sign(boot_ci_bonf)) > 0, 1, 0)

    results <- data.frame(obs_estimate=obs_est, unbiased_estimate=2 * obs_est - mean(boot_est), boot_ci, boot_ci_bonf, signif_ci=signif_ci, signif_ci.adj=signif_ci_adj)
  
    return(results)
  })
  
  pval_wild_boot <- data.frame(cluster=j, type=c(ps_no_intercept, ps_intercept), Reduce(rbind, wild_boot))

  return(pval_wild_boot)
})

wild_boot_type <- Reduce(rbind, list_wild_boot_type) %>% arrange(-signif_ci)

list_wild_boot_timepoint <- lapply(cluster_names, function(j) {
  cluster_data_perc <- subset(data_perc, cluster == j)
  cluster_data_perc$timepoint <- relevel(factor(cluster_data_perc$timepoint), ref=group_col_baseline)

  fit <- lm(value ~ type + timepoint, data=cluster_data_perc)
  
  timepoint_no_intercept <- na.omit(str_extract(names(coef(fit)), 'timepoint.*'))
  timepoint_intercept <- setdiff(timepoint_names, gsub('timepoint', '', timepoint_no_intercept))

  wild_boot <- lapply(c(timepoint_no_intercept, '(Intercept)'), function(i) {
    if (i != '(Intercept)') {
      cat('running wild bootstrap for', j, gsub('timepoint', '', i), '...\n')

      suppressWarnings({
          fit_wild_boot <- wild.boot(value ~ type + timepoint, data=cluster_data_perc, B=B)
        })

        obs_est <- coef(fit)[i]
        boot_est <- fit_wild_boot$bootEstParam[, i]

        limits <- c(alpha_level / k, 1 - alpha_level / k)
        boot_ci <- quantile(boot_est, probs=limits) %>% matrix(nrow=1)
        boot_ci[1] <- ifelse(test_kind == 'lower-tailed', -Inf, boot_ci[1])
        boot_ci[2] <- ifelse(test_kind == 'upper-tailed', Inf, boot_ci[2])
        colnames(boot_ci) %<>% paste0('ci_limit_', limits)
        colnames(boot_ci)[1] <- ifelse(test_kind == 'lower-tailed', 'ci_limit_minusInf', colnames(boot_ci)[1])
        colnames(boot_ci)[2] <- ifelse(test_kind == 'upper-tailed', 'ci_limit_Inf', colnames(boot_ci)[2])
        signif_ci <- ifelse(prod(sign(boot_ci)) > 0, 1, 0)

        adj <- length(timepoint_no_intercept) + 1
        limits_bonf <- c((alpha_level / adj) / k, 1 - (alpha_level / adj) / k)
        boot_ci_bonf <- quantile(boot_est, probs=limits_bonf) %>% matrix(nrow=1)
        boot_ci_bonf[1] <- ifelse(test_kind == 'lower-tailed', -Inf, boot_ci_bonf[1])
        boot_ci_bonf[2] <- ifelse(test_kind == 'upper-tailed', Inf, boot_ci_bonf[2])
        colnames(boot_ci_bonf) %<>% paste0('ci.adj_limit_', limits_bonf)
        colnames(boot_ci_bonf)[1] <- ifelse(test_kind == 'lower-tailed', 'ci.adj_limit_minusInf', colnames(boot_ci_bonf)[1])
        colnames(boot_ci_bonf)[2] <- ifelse(test_kind == 'upper-tailed', 'ci.adj_limit_Inf', colnames(boot_ci_bonf)[2])
        signif_ci_adj <- ifelse(prod(sign(boot_ci_bonf)) > 0, 1, 0)

        results <- data.frame(obs_estimate=obs_est, unbiased_estimate=2 * obs_est - mean(boot_est), boot_ci, boot_ci_bonf, signif_ci=signif_ci, signif_ci.adj=signif_ci_adj)
      } else {
        cat('appending values for', j, timepoint_intercept, '...\n')

        limits <- c(alpha_level / k, 1 - alpha_level / k)
        boot_ci <- matrix(c(0, 0), nrow=1)
        colnames(boot_ci) %<>% paste0('ci_limit_', limits)
        colnames(boot_ci)[1] <- ifelse(test_kind == 'lower-tailed', 'ci_limit_minusInf', colnames(boot_ci)[1])
        colnames(boot_ci)[2] <- ifelse(test_kind == 'upper-tailed', 'ci_limit_Inf', colnames(boot_ci)[2])

        adj <- length(timepoint_no_intercept) + 1
        limits_bonf <- c((alpha_level / adj) / k, 1 - (alpha_level / adj) / k)
        boot_ci_bonf <- matrix(c(0, 0), nrow=1)
        colnames(boot_ci_bonf) %<>% paste0('ci.adj_limit_', limits_bonf)
        colnames(boot_ci_bonf)[1] <- ifelse(test_kind == 'lower-tailed', 'ci.adj_limit_minusInf', colnames(boot_ci_bonf)[1])
        colnames(boot_ci_bonf)[2] <- ifelse(test_kind == 'upper-tailed', 'ci.adj_limit_Inf', colnames(boot_ci_bonf)[2])

        results <- data.frame(obs_estimate=0, unbiased_estimate=0, boot_ci, boot_ci_bonf, signif_ci=0, signif_ci.adj=0)
      }

      return(results)
    })
  
  pval_wild_boot <- data.frame(cluster=j, timepoint=gsub('timepoint', '', c(timepoint_no_intercept, timepoint_intercept)), Reduce(rbind, wild_boot))

  return(pval_wild_boot)
})

wild_boot_timepoint <- Reduce(rbind, list_wild_boot_timepoint) %>% arrange(-signif_ci)

# DATA VISUALIZATION -----------------------------------------------------------

list_boxplot_type <- lapply(cluster_names, function(j) {
  cluster_data_perc <- subset(data_perc, cluster == j)

  mean_type <- ddply(cluster_data_perc, ~ cluster + type, summarise, m=mean(value))
  
  if (!is.null(ordered_type)) {
    mean_type$type %<>% factor(levels=ordered_type, ordered=TRUE)
    cluster_data_perc$type %<>% factor(levels=ordered_type, ordered=TRUE)
  }

  p_boxplot_type <- ggplot(mean_type) +
    geom_boxplot(inherit.aes=FALSE, data=cluster_data_perc, aes(x=type, y=value), alpha=0.2) +
    geom_point(aes(x=type, y=m), alpha=0.8, color='red', size=2) +
    geom_path(aes(x=as.numeric(as.factor(type)), y=m), alpha=0.2, color='red', linewidth=1) +
    geom_hline(aes(yintercept=mean(m)), linetype='dashed') +
    facet_wrap(~ cluster, scales='free_y') +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.6)) +
    labs(y='Relative frequency', x='') +
    scale_y_continuous(labels = scales::percent)

  return(p_boxplot_type)
})
names(list_boxplot_type) <- cluster_names

list_boxplot_timepoint <- lapply(cluster_names, function(j) {
  cluster_data_perc <- subset(data_perc, cluster == j)

  mean_timepoint <- ddply(cluster_data_perc, ~ cluster + timepoint, summarise, m=mean(value))
  
  if (!is.null(ordered_timepoint)) {
    mean_timepoint$timepoint %<>% factor(levels=ordered_timepoint, ordered=TRUE)
    cluster_data_perc$timepoint %<>% factor(levels=ordered_timepoint, ordered=TRUE)
  }
  
  p_boxplot_timepoint <- ggplot(mean_timepoint) +
    geom_boxplot(inherit.aes=FALSE, data=cluster_data_perc, aes(x=timepoint, y=value), alpha=0.2) +
    geom_point(aes(x=timepoint, y=m), alpha=0.8, color='red', size=2) +
    geom_path(aes(x=as.numeric(as.factor(timepoint)), y=m), alpha=0.2, color='red', linewidth=2) +
    geom_hline(aes(yintercept=mean(m)), linetype='dashed') +
    facet_wrap(~ cluster, scales='free_y') +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.6)) +
    labs(y='Relative frequency', x='') +
    scale_y_continuous(labels = scales::percent)

  return(p_boxplot_timepoint)
})
names(list_boxplot_timepoint) <- cluster_names

if (!is.null(ordered_cluster)) wild_boot_type$cluster %<>% factor(levels=ordered_cluster, ordered=TRUE)

if (!is.null(ordered_type)) wild_boot_type$type %<>% factor(levels=ordered_type, ordered=TRUE)
if (p_kind == 'unadjusted') {
  wild_boot_type$signif_res <- wild_boot_type$signif_ci
} else {
  wild_boot_type$signif_res <- wild_boot_type$signif_ci.adj
}

p_heatmap_type <- ggplot(wild_boot_type, aes(x = cluster, y = type)) +
  geom_tile(aes(fill = unbiased_estimate), color = "black") +
  geom_point(aes(shape = ifelse(signif_res == 1, "y", "n")), size = 3,
             color = "black", show.legend = FALSE) +
  scale_size(range = c(1, 5)) +
  scale_fill_gradient2(low = "tomato3", mid = "white", high = "forestgreen") +
  scale_shape_manual(values = c(y = 8, n = NA)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "", size = "p.adj", fill = "Diff. from grand-mean:")

if (!is.null(ordered_timepoint)) wild_boot_timepoint$timepoint %<>% factor(levels=ordered_timepoint, ordered=TRUE)
if (p_kind == 'unadjusted') {
  wild_boot_timepoint$signif_res <- wild_boot_timepoint$signif_ci
} else {
  wild_boot_timepoint$signif_res <- wild_boot_timepoint$signif_ci.adj
}

p_heatmap_timepoint <- ggplot(wild_boot_timepoint, aes(x = cluster, y = timepoint)) +
  geom_tile(aes(fill = unbiased_estimate), color = "black") +
  geom_point(aes(shape = ifelse(signif_res == 1, "y", "n")), size = 3,
             color = "black", show.legend = FALSE) +
  scale_size(range = c(1, 5)) +
  scale_fill_gradient2(low = "tomato3", mid = "white", high = "forestgreen") +
  scale_shape_manual(values = c(y = 8, n = NA)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "", size = "p.adj", fill = paste0("Diff. from ", group_col_baseline, ":"))

# SAVE RESULTS -----------------------------------------------------------------

sys_time <- substr(gsub(':', '-', gsub(' ', '_', Sys.time())), 1, 19)
if (!dir.exists(paste0(path, '/outputs/'))) dir.create(paste0(path, '/outputs/'))
base::dir.create(base::paste0(path, "/outputs/", file_name, '_', sys_time, '_results/'))
dir.create(paste0(path, '/outputs/', file_name, '_', sys_time, '_results', '/boxplots'))
dir.create(paste0(path, '/outputs/', file_name, '_', sys_time, '_results', '/', test_kind))

lapply(seq_len(length(list_boxplot_type)), function(i) {
  fig <- list_boxplot_type[i]
  clust <- names(fig)
  ggsave(paste0(path, "/outputs/", file_name, '_', sys_time, "_results", "/boxplots/type_boxplot_mean_grand-mean_", clust, ".pdf"),
         fig[[clust]], dpi = 125, width = 45, height = 20, units = "cm")
})

lapply(seq_len(length(list_boxplot_timepoint)), function(i) {
  fig <- list_boxplot_timepoint[i]
  clust <- names(fig)
  ggsave(paste0(path, "/outputs/", file_name, '_', sys_time, "_results", "/boxplots/", group_col, "_boxplot_mean_grand-mean_", clust, ".pdf"),
         fig[[clust]], dpi = 125, width = 45, height = 20, units = "cm")
})

ggsave(paste0(path, "/outputs/", file_name, '_', sys_time, "_results", '/', test_kind, "/type_heatmap-all-clusters.pdf"),
       p_heatmap_type, dpi = 125, width = 20, height = 20, units = "cm")

ggsave(paste0(path, "/outputs/", file_name, '_', sys_time, "_results", '/', test_kind, "/", group_col, "_heatmap-all-clusters.pdf"),
       p_heatmap_timepoint, dpi = 125, width = 20, height = 20, units = "cm")

write_xlsx(list(type=dplyr::select(wild_boot_type, -signif_res),
                group = dplyr::select(wild_boot_timepoint, -signif_res)),
                path = paste0(path, "/outputs/", file_name, '_', sys_time, "_results", '/', test_kind, "/estimates-and-pvalues.xlsx"))

write_xlsx(list(data=data_perc, ps_kept_excluded=data_changes_outline),
           path = paste0(path, "/outputs/", file_name, '_', sys_time, "_results", "/dataset.xlsx"))

# OBJECTS HOLDING FINAL RESULTS ------------------------------------------------

data_perc # data

list_boxplot_type # named list with several boxplots for comparisons between PSs
list_boxplot_timepoint # named list with several boxplots for comparisons between timepoints

p_heatmap_type # heatmap for the difference from the grand-mean. comparisons between PSs
p_heatmap_timepoint # heatmap for the difference from the grand-mean. comparisons between timepoints

wild_boot_type # table of results for the wild bootstrap (estimates and p-values). comparisons between PSs
wild_boot_timepoint # table of results for the wild bootstrap (estimates and p-values). comparisons between timepoints

