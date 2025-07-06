library(dplyr)
library(purrr)
library(ggplot2)
library(ggnewscale)
library(gridExtra)

# Edge ratio
E <- seq(0, 1, length.out = 500)

# Global parameters
N_H <- 1
N_M0 <- 1
N_W0 <- 5
k_W <- 3
Q_H <- 0.3
Q_M0 <- 0.6
Q_W <- 0.1
f <- 1/3
c_M <- 0.27
D_max <- 1
E_p <- 0.4
sigma_mosq <- 0.15

# Mosquito density function
D_V <- D_max * exp(-((E - E_p)^2) / (2 * sigma_mosq^2))

# Baseline
N_M <- rep(N_M0, length(E))
N_W <- N_W0 * exp(-k_W * E)
total_bites <- N_H * Q_H + N_M * Q_M0 + N_W * Q_W
bite_H <- (N_H * Q_H) / total_bites
bite_M <- (N_M * Q_M0) / total_bites
C_E <- c_M * bite_M
H_E <- f * bite_H
R_baseline <- D_V * N_M * C_E * H_E
R_baseline_rel <- R_baseline / max(R_baseline)
baseline_df <- data.frame(EdgeRatio = E, Risk = R_baseline_rel, Param = "Baseline", Scenario = "Baseline")

# Scenario 1: Vary k_M
k_M_vals <- c(1, 2, 3)
s1 <- map_dfr(k_M_vals, function(k_M) {
  N_M1 <- N_M0 * exp(-k_M * E)
  total_bites <- N_H * Q_H + N_M1 * Q_M0 + N_W * Q_W
  bite_H <- (N_H * Q_H) / total_bites
  bite_M <- (N_M1 * Q_M0) / total_bites
  C_E <- c_M * bite_M
  H_E <- f * bite_H
  Risk <- D_V * N_M1 * C_E * H_E / max(R_baseline)
  data.frame(EdgeRatio = E, Risk = Risk, Param = paste0("k_M = ", k_M), Scenario = "Scenario 1")
})

# Scenario 2: Vary biodiversity collapse threshold
threshold_vals <- c(0.4, 0.5, 0.6)
s2 <- map_dfr(threshold_vals, function(thresh) {
  N_W2 <- ifelse(E < thresh, N_W0 * (1 - E), 0.001)
  total_bites <- N_H * Q_H + N_M * Q_M0 + N_W2 * Q_W
  bite_H <- (N_H * Q_H) / total_bites
  bite_M <- (N_M * Q_M0) / total_bites
  C_E <- c_M * bite_M
  H_E <- f * bite_H
  Risk <- D_V * N_M * C_E * H_E / max(R_baseline)
  data.frame(EdgeRatio = E, Risk = Risk, Param = paste0("Threshold = ", thresh), Scenario = "Scenario 2")
})

# Scenario 3: Vary human density
N_H_max_vals <- c(0.1, 2, 5, 20)
s3 <- map_dfr(N_H_max_vals, function(N_H_max) {
  N_H3 <- 1 + (N_H_max - 1) * E  # Linear increase from 1 to N_H_max
  total_bites <- N_H3 * Q_H + N_M * Q_M0 + N_W * Q_W
  bite_H <- (N_H3 * Q_H) / total_bites
  bite_M <- (N_M * Q_M0) / total_bites
  C_E <- c_M * bite_M
  H_E <- f * bite_H
  Risk <- D_V * N_M * C_E * H_E / max(R_baseline)
  data.frame(EdgeRatio = E, Risk = Risk, Param = paste0("N_H_max = ", N_H_max), Scenario = "Scenario 3")
})

# Combine and add baseline for all panels
plot_data <- bind_rows(s1, s2, s3)
plot_data <- bind_rows(plot_data)

################################################################################

plot_data_expanded <- plot_data %>%
  mutate(
    Param = case_when(
      Scenario == "Scenario 1" ~ paste0("k[M]==", gsub("k_M = ", "", Param)),
      Scenario == "Scenario 2" ~ paste0("E[thres]==", gsub("Threshold = ", "", Param)),
      Scenario == "Scenario 3" ~ paste0("N[H]^max==", gsub("N_H_max = ", "", Param)),
      TRUE ~ Param
    )
  ) %>%
  bind_rows(
    baseline_df %>% mutate(Scenario = "Scenario 1", Param = "Baseline"),
    baseline_df %>% mutate(Scenario = "Scenario 2", Param = "Baseline"),
    baseline_df %>% mutate(Scenario = "Scenario 3", Param = "Baseline")
  )

plot_data_expanded <- plot_data_expanded %>%
  mutate(Param = factor(Param, 
                        levels = c("Baseline", "k[M]==1", "k[M]==2", "k[M]==3",  # Scenario 1 order
                                   "E[thres]==0.4", "E[thres]==0.5", "E[thres]==0.6",  # Scenario 2
                                   "N[H]^max==0.1", "N[H]^max==2", "N[H]^max==5", "N[H]^max==20"),  # Scenario 3
                        ordered = TRUE)
  )

# Create color gradients for each scenario
scenario_colors <- list(
  "Scenario 1" = scale_color_manual(
    values = c("#1b9e77", "#d95f02", "#7570b3"),
    labels = function(x) gsub(".*=", "", x), 
    name = expression(italic(k[M]))
  ),
  "Scenario 2" = scale_color_manual(
    values = c("#e7298a", "#66a61e", "#e6ab02"),
    labels = function(x) gsub(".*=", "", x),
    name = expression(italic(Epsilon[thres]))
  ),
  "Scenario 3" = scale_color_manual(
    values = c('gray','skyblue','brown','green'),
    labels = function(x) gsub(".*=", "", x),
    name = expression(italic(N[H]^max))
  )
)

# Create the plot
a=ggplot() + ggtitle('(a) Scenario 1')+
  # Scenario 1
  geom_line(
    data = filter(plot_data_expanded, Scenario == "Scenario 1" & Param != "Baseline"),
    aes(x = EdgeRatio, y = Risk, color = Param)
  )  +
  scenario_colors[["Scenario 1"]] +
  geom_line(data = filter(plot_data_expanded, Param == "Baseline"),aes(x = EdgeRatio, y = Risk),color = "black",linetype = "dashed",alpha = 0.7) +
  labs(x = "",y = expression("" * italic(R(Epsilon)) * "")) +
  theme_minimal() +
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"),
        axis.text.x = element_blank(),
        legend.position = "right",
        legend.text.align = 0,
        legend.box.margin = margin(0,20,0,0),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(size = 8)
  )  + ylim(0,1.4)

# Scenario 2
b=ggplot() +  new_scale_color() +  ggtitle('(b) Scenario 2')+
  geom_line(
    data = filter(plot_data_expanded, Scenario == "Scenario 2" & Param != "Baseline"),
    aes(x = EdgeRatio, y = Risk, color = Param)) +
  scenario_colors[["Scenario 2"]] +
  geom_line(data = filter(plot_data_expanded, Param == "Baseline"),aes(x = EdgeRatio, y = Risk),color = "black",linetype = "dashed",alpha = 0.7) +
  labs(x = "",y = expression("" * italic(R(Epsilon)) * "")) +
  theme_minimal() +
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"),
        axis.text.x = element_blank(),
        legend.position = "right",
        legend.text.align = 0,
        legend.box.margin = margin(0,0,0,0),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(size = 8)
  )  + ylim(0,1.4)

# Scenario 3 
c=ggplot()+  new_scale_color() +  ggtitle('(c) Scenario 3')+
  geom_line(
    data = filter(plot_data_expanded, Scenario == "Scenario 3" & Param != "Baseline"),
    aes(x = EdgeRatio, y = Risk, color = Param)) +
  scenario_colors[["Scenario 3"]] +
  geom_line(data = filter(plot_data_expanded, Param == "Baseline"),aes(x = EdgeRatio, y = Risk),color = "black",linetype = "dashed",alpha = 0.7) +
  labs(x = expression("Edge Ratio (" * italic(Epsilon) * ")"),y = expression("" * italic(R(Epsilon)) * "")) +
  theme_minimal() +
  theme(plot.margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"),
        legend.position = "right",
        legend.text.align = 0,
        legend.box.margin = margin(0,0,0,0),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(size = 8)
  ) + ylim(0,1.4)

f3=grid.arrange(a,b,c)

ggsave("Fig3.png", 
       f3, width = 85, height = 130, units = c("mm"), dpi = 300)
