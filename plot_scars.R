library(argparse)
library(karyoploteR)

# ------------------------------ # 
parser <- ArgumentParser(description="Visualise HRD scars")
parser$add_argument('input')
args <- parser$parse_args()


# parse HRD output
hrd_segments <- read.csv(args$input, sep="\t")
loh_seg <- hrd_segments[hrd_segments$HRD_tag == 'LOH', 
                        c("Chromosome", "Start_position", "End_position")]

tai_seg <- hrd_segments[hrd_segments$HRD_tag == 'TAI', 
                        c("Chromosome", "Start_position", "End_position")]

lst_seg <- hrd_segments[hrd_segments$HRD_tag == 'LST', 
                        c("Chromosome", "Start_position", "End_position")]

colnames(loh_seg) <- c("chr", "start", "end")
colnames(tai_seg) <- c("chr", "start", "end")
colnames(lst_seg) <- c("chr", "start", "end")


# plot HRD scars
ofname <- paste(strsplit(args$input, "[.]")[[1]][1], ".hrd.png", sep="")
png(ofname, width = 1400, height = 800)

kp <- plotKaryotype(genome = "hg19", chromosomes="autosomal")
kpPlotRegions(kp, data = loh_seg, col="#B2182B", r0=0, r1=0.25)
kpPlotRegions(kp, data = tai_seg, col="#2166AC", r0=0.3, r1=0.55)
kpPlotRegions(kp, data = lst_seg, col="#9933CC", r0=0.6, r1=0.85)

legend("bottomright", 
        legend=c("LOH", "TAI", "LST"), 
        col=c("#B2182B", "#2166AC", "#9933CC"),
        pch=c(15, 15, 15))

dev.off()
