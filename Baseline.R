# Load libraries once
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# 1. Define Epsilon (edge ratio) instead of epsilon
Epsilon <- seq(0, 1, length.out = 500)

# 2. Parameters ------------------------------------------------------------
N_H  <- 1      # Normalized human density
N_M  <- 1      # Constant macaque density
N_W0 <- 5      # Baseline wildlife density at Epsilon = 0
k_w  <- 3      # Wildlife‐loss rate

# Mosquito biting preferences
Q_H <- 0.3
Q_M <- 0.6
Q_W <- 0.1

# Biting frequency and macaque competence
f   <- 1/3
c_M <- 0.27

# Mosquito‐density Gaussian parameters
D_max     <- 1
Epsilon_p     <- 0.4
sigma     <- 0.15

# 3. Compute each component as a function of Epsilon --------------------------

# Wildlife density: N_W(Epsilon) = N_W0 * exp(-k_w * Epsilon)
N_W <- N_W0 * exp(-k_w * Epsilon)

# Mosquito density: D_V(Epsilon) = D_max * exp(-((Epsilon - Epsilon_p)^2)/(2 * sigma^2))
D_V <- D_max * exp(-((Epsilon - Epsilon_p)^2)/(2 * sigma^2))

# Biting‐allocation denominators
denom <- (N_H * Q_H) + (N_M * Q_M) + (N_W * Q_W)

# Proportion of bites on each host
bite_H <- (N_H * Q_H) / denom
bite_M <- (N_M * Q_M) / denom
bite_W <- (N_W * Q_W) / denom

# Community competence: C(Epsilon) = c_M * (bites on macaques)
C_Epsilon <- c_M * bite_M

# Human contact rate: H(Epsilon) = f * (bites on humans)
H_Epsilon <- f * bite_H

# Unnormalized spillover risk: R(Epsilon) = D_V * N_M * C(Epsilon) * H(Epsilon)
R_Epsilon <- D_V * N_M * C_Epsilon * H_Epsilon

# Relative‐risk scaling
R_rel <- R_Epsilon / max(R_Epsilon)

# 4. Combine into a data frame (with column named 'Epsilon')
results <- data.frame(
  Epsilon                 = Epsilon,
  WildlifeDensity     = N_W,
  MosquitoDensity     = D_V,
  CommunityCompetence = C_Epsilon,
  HumanContactRate    = H_Epsilon,
  SpilloverRisk       = R_rel
)

# 5. Reshape for faceting and prepend (a)–(e) to each strip label -----------
results_long <- results %>%
  pivot_longer(cols = -Epsilon, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric =
           recode(Metric,
                  "MosquitoDensity"      = "(a)~D[v](Epsilon)",
                  "WildlifeDensity"      = "(b)~N[W](Epsilon)",
                  "CommunityCompetence"  = "(c)~C(Epsilon)",
                  "HumanContactRate"     = "(d)~H(Epsilon)",
                  "SpilloverRisk"        = "(e)~R(Epsilon)"
           )
  ) %>%
  mutate(Metric = factor(Metric, levels = c(
    "(a)~D[v](Epsilon)",
    "(b)~N[W](Epsilon)",
    "(c)~C(Epsilon)",
    "(d)~H(Epsilon)",
    "(e)~R(Epsilon)"
  )))

# 6. Plot with facet labels parsed (so Greek Epsilon is rendered, and each strip shows (A)–(E)) -------------
f2 <- ggplot(results_long, aes(x = Epsilon, y = Value, group = Metric)) +
  geom_line(aes(color = Metric)) +
  facet_grid(
    rows = vars(Metric),
    scales = "free_y",
    switch = "y",
    labeller = labeller(Metric = label_parsed)
  ) +
  labs(
    x = expression("Edge Ratio (" * italic(Epsilon) * ")"),
    y = NULL
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 6)) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal(base_size = 12) +
  theme(
    axis.line.x       = element_line(color = "black", size = 0.3),
    axis.title.y = element_text(face = "italic"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11, face = "italic"),
    strip.background  = element_blank(),
    panel.spacing     = unit(0.8, "lines"),
    legend.position   = "none",
    axis.line         = element_line(color = "grey60"),
    axis.text         = element_text(color = "grey40"),
    axis.title.x      = element_text(color = "grey30", margin = margin(t = 8)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("Fig2.png", 
       f2, width = 85, height = 120, units = c("mm"), dpi = 300)
