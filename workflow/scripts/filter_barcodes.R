# Written by P.A. Gagnaire and modified by Alexis Simon.
# Filters the 10X barcodes on their coverage, 
# removing low coverage ones (before the first minimum of the function)
# and the high coverage one (after 10*MAX)

library(ggplot2)
library(stats)
library(scales)

barcode_input <- snakemake@input[[1]]
barcode_output <- snakemake@output[["barcodes"]]
fig_out <- snakemake@output[["figure"]]

X <- read.table(barcode_input, header = F)

H <- hist(X$V2, breaks = seq(0, max(X$V2), 1), plot = F)
D <- data.frame(cbind(H$breaks[2:(max(X$V2)+1)], H$counts))
NZ <- D[!(D$X2 %in% 0), ]
TOT <- sum((NZ$X1)*(NZ$X2))
count_empirical <- splinefun(x = NZ$X1, y = NZ$X2)
MIN <- floor(optimize(function(x) count_empirical(x), c(1, 50))$minimum)
MAX <- round((optimize(function(x) count_empirical(x), c(MIN, 500), maximum = T))$maximum)
LOW_REM <- sum((NZ[NZ$X1 < MIN, ]$X1)*(NZ[NZ$X1 < MIN, ]$X2))
UP_REM <- sum((NZ[NZ$X1 > 10*MAX, ]$X1)*(NZ[NZ$X1 > 10*MAX, ]$X2))
NKEPT <- TOT - LOW_REM - UP_REM
PKEPT <- NKEPT/TOT
X_out <- X[X$V2 >= MIN & X$V2 <= 10*MAX, ]

P <- ggplot(NZ, aes(x = X1, y = X2)) + xlab("Read pairs in bacode") + ylab("Count") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) +
  annotate("rect", xmin = MIN, xmax = 10*MAX, ymin = 0, ymax = +Inf, fill = 'green', alpha = 0.2) +
  annotate("rect", xmin = 0, xmax = MIN, ymin = 0, ymax = +Inf, fill = 'red', alpha = 0.2) + 
  annotate("rect", xmin = 10*MAX, xmax = +Inf, ymin = 0, ymax = +Inf, fill='red', alpha = 0.2) + 
  annotate("text", x = MAX, y = 5e4, label = "Maximum") + 
  annotate("text", x = MAX, y = 3e4, label = MAX) +
  annotate("text", x = MAX, y = 1.7e5, label = "Mean") + 
  annotate("text", x = MAX, y = 1e5, label = round(mean(X_hi$V2), 1)) +
  annotate("text", x = 3, y = 100, label = "Removed low") + 
  annotate("text", x = MAX, y = 100, label = "Retained") + 
  annotate("text", x = 2000,y = 100, label = "Removed high") +
  annotate("text", x = 3, y = 55, label = LOW_REM) + 
  annotate("text", x = MAX, y = 55, label = NKEPT) + 
  annotate("text", x = 2000, y = 55, label = UP_REM) +
  annotate("text", x = 3, y = 30, label = round(LOW_REM/TOT, 3)) + 
  annotate("text", x = MAX, y = 30, label = round(NKEPT/TOT, 3)) + 
  annotate("text", x = 2000, y = 30, label = round(UP_REM/TOT, 3)) +
  geom_point() + geom_line(color = "blue") + theme_bw() 

ggsave(fig_out, P, width = 8, height = 4)
writeLines(X_out$V1, barcode_output)

