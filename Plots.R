library(tidyverse)
library(ggthemes)
pal = tableau_color_pal()(10) #color palette
pal[c(7,1)]
54
unsup.df <- read.csv("unsup_full_df.csv")
unsup.half.df <- read.csv("unsup_half_df.csv")
unsup.hvg.df <- read.csv("unsup_hvg_df.csv")

colnames(unsup.half.df) <- colnames(unsup.df)
colnames(unsup.hvg.df) <- colnames(unsup.df)

## function to plot BIC
plot_bic <- function(df) {
        bic_df <- df[,c(1, 8:9)]
        bic_long_df <- bic_df %>% gather(cohort_type, bic, no_cohort_bic:if_cohort_bic, factor_key = TRUE) %>% 
                mutate(cohort = ifelse(cohort_type == "no_cohort_bic", "no cohort", "IF cohort"))
        
        highlight_df <- bic_long_df %>% 
                filter(n_starts == 10)
        
        bic_plot <- ggplot(bic_long_df, aes(x = as.factor(n_starts), y = bic, color = cohort)) +
                geom_point() +
                scale_color_manual(values = pal[c(7,1)]) +
                geom_smooth(aes(group = cohort), method = "loess", size = 2.5) +
                geom_point(data = highlight_df, 
                           aes(x = as.factor(n_starts), y = bic), 
                           color = 'black',
                           pch = 9,
                           size = 3.5) +
                ylab("BIC") +
                xlab("Number of iterations") +
                theme_classic()
        
        return(bic_plot)
        }

## function to plot computation time
plot_time <- function(df) {
        time_df <- df[,c(1, 10:11)]
        time_long_df <- time_df %>% 
                gather(cohort_type, time, noco_time:ifco_time, factor_key = TRUE) %>% 
                mutate(cohort = ifelse(cohort_type == "noco_time", "no cohort", "IF cohort"))
       
         highlight_df <- time_long_df %>% 
                filter(n_starts == 10)
        
        time_plot <- ggplot(time_long_df, aes(x = as.factor(n_starts), y = time, color = cohort)) +
                geom_point()+
                geom_smooth(aes(group = cohort), size = 2.5) +
                geom_point(data = highlight_df, 
                           aes(x = as.factor(n_starts), y = time), 
                           color = 'black',
                           pch = 9,
                           size = 3.5) +
                scale_color_manual(values = pal[c(7,1)]) +
                ylab("Computation") +
                xlab("Number of iterations") +
                theme_classic()
        
        
        return(time_plot)
        }

## function to plot ARI
plot_ari <- function(df) {
        ari_df <- df[,c(1:7)]
        ari_long_df <- (ari_df %>% 
                                gather(key = "cohort_type", value = "ari", no_cohort_ari, if_cohort_ari) %>% 
                                gather(key = "ci_lower", value = "ari_ci_lower", noco_ari_ci_lower, ifco_ari_ci_lower) %>%
                                gather(key = "ci_upper", value = "ari_ci_upper", noco_ari_ci_upper, ifco_ari_ci_upper)) %>% 
                mutate(cohort = ifelse(cohort_type == "no_cohort_ari", "no cohort", "IF cohort")) %>% 
                select(-c(ci_lower, ci_upper))
        
        highlight_df <- ari_long_df %>% 
                filter(n_starts == 10)
        
        ari_plot <- ggplot(ari_long_df, aes(x = as.factor(n_starts), y = ari, color = cohort)) +
                geom_point()+
                geom_smooth(aes(group = cohort), size = 2.5) +
                geom_point(data = highlight_df, 
                           aes(x = as.factor(n_starts), y = ari), 
                           color = 'black',
                           pch = 9,
                           size = 3.5) +
                scale_color_manual(values = pal[c(7,1)]) +
                ylab("ARI") +
                xlab("Number of iterations") +
                theme_classic()
        
        return(ari_plot)
        }

# create a list of data frames and a list of function names
df_list <- list(unsup.df, unsup.half.df, unsup.hvg.df)
func_list <- list("plot_ari", "plot_bic", "plot_time")

# create an empty list to store the output plots
output_plots <- list()

# iterate over the data frames and apply the plotting functions
for (i in seq_along(df_list)) {
        for (j in seq_along(func_list)) {
                # get the current data frame and function name
                df <- df_list[[i]]
                func <- func_list[[j]]
                
                # generate the plot and save it with a unique name
                plot <- do.call(func, list(df))
                filename <- paste0("plot_nstart", i, "_", func, ".png")
                ggsave(filename, plot, device = "png")
                
                # add the plot to the output_plots list
                output_plots[[length(output_plots) + 1]] <- plot
        }
}





