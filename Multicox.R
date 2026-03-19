cox_model_3 <- coxph(Surv(OS.time, OS) ~  Age + Stage + PCDI + IRGI + ICDI, data = data)

# Extract and summarize Cox model
multivar_all_3 <- summary(cox_model_3)

# Extract components
HR <- multivar_all_3$coefficients[, "exp(coef)"]
CI_lower <- multivar_all_3$conf.int[, "lower .95"]
CI_upper <- multivar_all_3$conf.int[, "upper .95"]
p_values <- multivar_all_3$coefficients[, "Pr(>|z|)"]

# Combine into dataframe
result <- data.frame(
  Variable = rownames(multivar_all_3$coefficients),
  HR = HR,
  CI_lower_95 = CI_lower,
  CI_upper_95 = CI_upper,
  p_value = p_values
)

# Add formatted HR (95% CI) column
result$`HR (95% CI for HR)` <- paste0(
  sprintf("%.2f", result$HR),
  " (",
  sprintf("%.2f", result$CI_lower_95),
  ", ",
  sprintf("%.2f", result$CI_upper_95),
  ")"
)

# Now generate final output table
dt_all_3 <- data.frame(
  Variables = c("Age", "Stage", "PCDI", "IRGI","ICDI" ),
  P_value = ifelse(result$p_value < 0.001, "<0.001", signif(result$p_value, 2)),
  HR_CI = result$`HR (95% CI for HR)`,
  est = result$HR,
  low = result$CI_lower_95,
  hi = result$CI_upper_95
)

print(dt_all_3)

# Adjust the column widt_all_3h with space.
dt_all_3$` ` <- paste(rep(" ", 20), collapse = " ")

# Create confidence interval column to display
dt_all_3$`HR (95% CI)` <- dt_all_3$HR_CI

dt_all_3$hi <- as.numeric(dt_all_3$hi)
dt_all_3$est <- as.numeric(dt_all_3$est)
dt_all_3$low <- as.numeric(dt_all_3$low)

library(extrafont)

# Import fonts into R
#font_import()

# Load the font
#loadfonts(device = "win")

# Define theme
tm <- forest_theme(base_size = 10,
                   base_family = "Times New Roman", # Set font to Times New Roman
                   refline_col = "red",
                   arrow_type = "closed",
                   panel_bg = c("white", "#1F3B73")
)

p <- forest(dt_all_3[,c(1:2, 7:8)],
            est = dt_all_3$est,
            lower = dt_all_3$low,
            upper = dt_all_3$hi,
            ci_column = 3,
            ref_line = 1,
            arrow_lab = c("Low Risk", "High Risk"),
            xlim = c(0,3 ),
            ticks_at = c(0.5, 1, 2,3),
            theme = tm)
# Print plot
plot(p)
library(ggplot2)
# Save the plot with 300 dpi
ggsave("forest_plot_multivariate_all_3.tiff", plot = p, dpi = 600, width = 8, height = 6)




