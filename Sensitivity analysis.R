library(sensobol)
library(ggplot2)
library(dplyr)

# Generate the sensitivity indices 
N <- 2000
params <- c("N_W0", "k_w", "c_M", "epsilon_p", "sigma")
#set.seed(123)

scale_params <- function(mat) {
  mat_scaled <- mat
  mat_scaled[, "N_W0"] <- qunif(mat[, "N_W0"], min = 1, max = 10)
  mat_scaled[, "k_w"] <- qunif(mat[, "k_w"], min = 1, max = 5)
  mat_scaled[, "c_M"] <- qunif(mat[, "c_M"], min = 0.1, max = 0.5)
  mat_scaled[, "epsilon_p"] <- qunif(mat[, "epsilon_p"], min = 0.2, max = 0.6)
  mat_scaled[, "sigma"] <- qunif(mat[, "sigma"], min = 0.05, max =0.4)
  return(mat_scaled)
}

mat_raw <- sobol_matrices(N = N, params = params)
mat <- scale_params(mat_raw)

spillover_model <- function(X) {
  apply(X, 1, function(params) {
    N_W0 <- params[1]
    k_w <- params[2]
    c_M <- params[3]
    epsilon_p <- params[4]
    sigma <- params[5]
    
    N_H <- 1; N_M <- 1; f <- 1/3; D_max <- 1
    Q_H <- 0.3; Q_M <- 0.6 ;Q_W <- 0.1        
    epsilon <- seq(0, 1, length.out = 500)
    N_W <- N_W0 * exp(-k_w * epsilon)
    D_V <- D_max * exp(-( (epsilon - epsilon_p)^2 ) / (2 * sigma^2))
    total_biting_pref <- (N_H * Q_H) + (N_M * Q_M) + (N_W * Q_W)
    bite_H <- (N_H * Q_H) / total_biting_pref
    bite_M <- (N_M * Q_M) / total_biting_pref
    C_E <- c_M * bite_M
    H_E <- f * bite_H
    R_E <- D_V * N_M * C_E * H_E
    max(R_E)
  })
}

Y <- spillover_model(mat)
R <- 100
ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = R)

# Create publication-quality plots
# Extract first-order and total-order indices
first_order <- ind$results[ind$results$sensitivity == "Si", ]
total_order <- ind$results[ind$results$sensitivity == "Ti", ]

# Combine first-order and total-order indices into one dataframe
all_indices <- rbind(
  cbind(first_order, Index_Type = "First-order (Si)"),
  cbind(total_order, Index_Type = "Total-order (STi)")
)

# Create parameter labels with Greek letters
param_labels <- c(
  "N_W0" = expression(N[W0]),
  "k_w" = expression(kappa[w]),
  "c_M" = expression(c[M]),
  "epsilon_p" = expression(Epsilon[p]),
  "sigma" = expression(sigma)
)

# Create single combined plot
combined_plot <- ggplot(all_indices, 
                        aes(x = reorder(parameters, -original), 
                            y = original, 
                            fill = Index_Type)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.8),
           width = 0.7,
           alpha = 0.8) +
  geom_errorbar(aes(ymin = low.ci, ymax = high.ci),
                position = position_dodge(width = 0.8),
                width = 0.25,
                color = "#333333") +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"),
                    labels=c("First-order (Si)"="First-order", "Total-order (STi)"="Total-order")) +
  scale_x_discrete(labels = param_labels) +  # Apply Greek letter labels
  labs(x = "",
       y = "Sobol Sensitivity Index",
       fill = "Index Type") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 10, face='italic'),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

# Display the plot
print(combined_plot)

ggsave("Fig4.png", 
       combined_plot, width = 174, height = 80, units = c("mm"), dpi = 300, bg = "white")
