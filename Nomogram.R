library(rms)
library(survival)

dt <- data.frame(
  Barcodes = data_clean$Barcode,
  Age = data_clean$Age,
  ICDI = data_clean$ICDI,
  Stage = data_clean$Stage,
  time = data_clean$time,
  status = data_clean$status
)
dd <- datadist(dt)
options(datadist = "dd")
f <- cph(Surv(time, status) ~ Age + Stage + ICDI, data = dt,
         x = TRUE, y = TRUE, surv = TRUE, time.inc = 365.25)
surv <- Survival(f)

nom <- nomogram(f, 
                fun = list( 
                  function(x) surv(365.25, x),     
                  function(x) surv(365.25 * 3, x),
                  function(x) surv(365.25 * 5, x)
                ),
                lp = FALSE,
                funlabel = c("1-year survival", "3-year survival", "5-year survival"),
                maxscale = 10,
                fun.at = pretty(c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5)))
plot(nom)

# Save with high-res and Times New Roman
tiff("nomogram_custom_long.tiff", width = 1900, height = 1500, res = 300, family = "Times New Roman")

# Plot nomogram
plot(nom, lplabel = "Total Points", 
     col.grid = gray(0.9),  # grey reference lines
     col.axis = "black", 
     col.label = "black",
     col = "black", 
     xfrac = 0.4,
     label.every = 1,
     cex.axis = 0.6,
     cex.var = 0.7)


dev.off()
