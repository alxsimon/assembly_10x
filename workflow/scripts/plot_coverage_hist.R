library(ggplot2)
library(ggforce)
library(glue)
theme_set(theme_minimal(base_size = 12))

defaultW <- getOption("warn")
options(warn = -1)

# get params
in_hist = snakemake@input[['hist']]
out_pdf = snakemake@output[['pdf']]
sample_name = paste0(snakemake@wildcards[['sample']], "_", snakemake@wildcards[['version']])

df = read.table(in_hist, col.names = c("bin", "count"))[-1,]
df$bin_int = as.numeric(df$bin)
df$bin_int[1001] = 1001

# gross mean
m = sum(df$bin_int*df$count)/sum(df$count)
# take care of too low coverage cases
if (m >= 5){
	max_bin = 2*m
} else {
	max_bin = 5
}

p = ggplot(df, aes(x = bin_int, y = count)) + 
	geom_bar(stat = 'identity') +
	scale_x_continuous(expand = c(0,0)) +
	scale_y_continuous(expand = c(0,0)) +
	labs(
		x = "Coverage",
		y = "Count",
		title = glue("Histogram of coverage for {sample_name} (last bin is >=1000)")
	) +
	theme(
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		plot.margin = margin(10, 25, 10, 25),
		axis.ticks.x = element_line(),
		axis.ticks.length.x = unit(5, "pt"),
		strip.background = element_rect(fill = "grey90", linetype = "blank")
	) +
	facet_zoom(xy = (bin_int > 0 & bin_int <= max_bin), zoom.size = 1, horizontal = F)

ggsave(out_pdf, p, width = 10, height = 8)

options(warn = defaultW)