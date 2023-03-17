library(tidyverse)
library(ggthemes)
pal = tableau_color_pal()(10) #color palette
pal[c(7,1)]

### Plot Funcions ###

#input nocohorts df first and ifcohorts df second

#this will work for any benchmark parameter tested given data
#is saved in similar format to csv files saved from 
# 'nphase1_benchmark.qmd'

plot_bic <- function(nocohorts, ifcohorts) {
  
  bic_df_noco <- nocohorts[,c(1, 3)] 
  bic_df_noco$cohort <- c("no cohorts")
  bic_df_ifco <- ifcohorts[,c(1, 3)]
  bic_df_ifco$cohort <- c("IF cohorts")
  
  bic_long_df <- rbind(bic_df_noco, bic_df_ifco)
  
  highlight_df <- bic_long_df %>% filter(n_phase1 == 10000)
  
  bic_plot <- ggplot(bic_long_df, aes(x = as.factor(n_phase1), y = bic, color = cohort)) +
    geom_point(size = 3.5) +
    scale_color_manual(values=pal[c(7,1)]) +
    geom_smooth(aes(group = cohort), method = "loess", se=TRUE, size=2.5) +
    geom_point(data = highlight_df, 
               aes(x = as.factor(n_phase1), y = bic), 
               color = 'black',
               pch = 9,
               size = 4) +
    ylab("BIC") +
    xlab("Phase 1 cell size") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(bic_plot)
}

#plot_bic(nphase1_full_gene_no_cohort, nphase1_full_gene_if_cohort)

plot_time <- function(nocohorts, ifcohorts) {
  
  df_noco <- nocohorts[,c(1, 4)] 
  df_noco$cohort <- c("no cohorts")
  df_ifco <- ifcohorts[,c(1, 4)]
  df_ifco$cohort <- c("IF cohorts")
  long_df <- rbind(df_noco, df_ifco)
  
  highlight_df <- long_df %>% filter(n_phase1 == 10000)
  
  time_plot <- ggplot(long_df, aes(x = as.factor(n_phase1), y = elapsed_time, color = cohort)) +
    geom_point(size = 3.5) +
    scale_color_manual(values=pal[c(7,1)]) +
    geom_smooth(aes(group = cohort), method = "loess", se=TRUE, size=2.5) +
    geom_point(data = highlight_df, 
               aes(x = as.factor(n_phase1), y = elapsed_time), 
               color = 'black',
               pch = 9,
               size = 4) +
    ylab("Computation time (minutes)") +
    xlab("Phase 1 cell size") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(time_plot)
}

#plot_time(nphase1_full_gene_no_cohort, nphase1_full_gene_if_cohort)

plot_ari <- function(nocohorts, ifcohorts) {
  
  df_noco <- nocohorts[,c(1, 2)] 
  df_noco$cohort <- c("no cohorts")
  df_ifco <- ifcohorts[,c(1, 2)]
  df_ifco$cohort <- c("IF cohorts")
  long_df <- rbind(df_noco, df_ifco)
  
  highlight_df <- long_df %>% filter(n_phase1 == 10000)
  
  ari_plot <- ggplot(long_df, aes(x = as.factor(n_phase1), y = ari, color = cohort)) +
    geom_point(size = 3.5) +
    scale_color_manual(values=pal[c(7,1)]) +
    geom_smooth(aes(group = cohort), method = "loess", se=TRUE, size=2.5) +
    geom_point(data = highlight_df, 
               aes(x = as.factor(n_phase1), y = ari), 
               color = 'black',
               pch = 9,
               size = 4) +
    ylab("ARI") +
    xlab("Phase 1 cell size") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(ari_plot)
}

#plot_ari(nphase1_full_gene_no_cohort, nphase1_full_gene_if_cohort)

# create a list of data frames and a list of function names
df_nocolist <- list(nphase1_full_gene_no_cohort) # insert no cohort df
df_ifcolist <- list(nphase1_full_gene_if_cohort) # insert if cohort df
func_list <- list("plot_ari", "plot_bic", "plot_time")

# create an empty list to store the output plots
output_plots <- list()

# iterate over the data frames and apply the plotting functions
for (i in seq_along(df_nocolist) & seq_along(df_ifcolist)) {
  for (j in seq_along(func_list)) {
    # get the current data frame and function name
    dfnoco <- df_nocolist[[i]]
    dfifco <- df_ifcolist[[i]]
    func <- func_list[[j]]
    
    # generate the plot and save it with a unique name
    plot <- do.call(func, list(dfnoco,dfifco))
    filename <- paste0("plot_nphase1", i, "FULL", func, ".png")
    ggsave(filename, plot, device = "png")
    
    # add the plot to the output_plots list
    output_plots[[length(output_plots) + 1]] <- plot
  }
}



