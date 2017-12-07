library(ggplot2)
library(pheatmap)

args <- commandArgs(trailingOnly=TRUE)

data_ref <- read.csv(args[1], sep=" ", header=FALSE, row.names=1)
data <- read.csv(args[2], sep=" ", header=FALSE, row.names=1)

my_name_ref <- args[3]
my_name <- args[4]

set.seed(123)

data_ref[1:5,1:5]
data_ref <- data_ref[order(data_ref$V2, decreasing=TRUE),]
data_ref[1:5,1:5]
data[1:5,1:5]
data <- data[rownames(data_ref),]
data[1:5,1:5]

data_ref_profile <- data_ref[2:ncol(data_ref)]
data_profile <- data[2:ncol(data)]

tiff(
	my_name_ref,
	width = 2*1000,
	height = 4*1000,
	units = "px",
	compression = "lzw",
	res = 300
)

my_breaks <- seq(0,2,length.out=20)
my_palette <- colorRampPalette(c("white","red"))(n = (length(my_breaks)-1))

data_ref_profile[data_ref_profile > 2] <- 2
data_profile[data_profile > 2] <- 2

data_ref_profile[1:5,1:5]
data_profile[1:5,1:5]
means.data_ref_profile <- colMeans(data_ref_profile)
means.data_profile <- colMeans(data_profile)

pheatmap(
	as.matrix(data_ref_profile),
	cluster_rows=FALSE,
	cluster_cols=FALSE,
	show_rownames=FALSE,
	show_colnames=FALSE,
	scale="none",
	breaks=my_breaks,
# 	annotation_row=row_annotation,
	main="ATAC Profile Around TSS\nLoading All Reads",
	border_color="none",
	color=my_palette
)
dev.off()

tiff(
	my_name,
	width = 2*1000,
	height = 4*1000,
	units = "px",
	compression = "lzw",
	res = 300
)

pheatmap(
	as.matrix(data_profile),
	cluster_rows=FALSE,
	cluster_cols=FALSE,
	show_rownames=FALSE,
	show_colnames=FALSE,
	scale="none",
	breaks=my_breaks,
# 	annotation_row=row_annotation,
	main="ATAC Profile around TSS\tLoading Overlapping Reads",
	border_color="none",
	color=my_palette
)
dev.off()
