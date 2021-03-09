library(argparse)
library(karyoploteR)
library(dplyr)


# ------------------------------ # 
parser <- ArgumentParser(description = "Visualise HRD scars")
parser$add_argument('segment', help = "The output HRD segment file")
parser$add_argument('seqz', help = "The binned seqz file")
args <- parser$parse_args()


# parse HRD output
hrd_segments <- read.csv(args$segment, sep="\t")
loh_seg <- hrd_segments[hrd_segments$HRD_tag == 'LOH', 
                        c("Chromosome", "Start_position", "End_position")]

tai_seg <- hrd_segments[hrd_segments$HRD_tag == 'TAI', 
                        c("Chromosome", "Start_position", "End_position")]

lst_seg <- hrd_segments[hrd_segments$HRD_tag == 'LST', 
                        c("Chromosome", "Start_position", "End_position")]

colnames(loh_seg) <- c("chr", "start", "end")
colnames(tai_seg) <- c("chr", "start", "end")
colnames(lst_seg) <- c("chr", "start", "end")


# get coverage and BAF data
seqz <- gzfile(args$seqz, "rt")
seqz <- read.table(seqz, header = T)

seqz <- seqz %>%
  dplyr::mutate(logR = log2(depth.ratio)) %>%
  dplyr::select(chromosome, position, logR, Af, Bf)


# ----------  plot  ----------- #
# ----------------------------- #
# ofname <- paste0(strsplit(args$segment, "[.]")[[1]][1], ".hrd.pdf")
# pdf(ofname)

for (chrom in seq(1,22)){
    chrom <- as.character(chrom)
    
    ofname <- paste0(strsplit(args$segment, "[.]")[[1]][1], "_chr", chrom, "_hrd.png")
    png(ofname, width = 1400, height = 800)
  
    plot.params <- getDefaultPlotParams(plot.type=2)
    plot.params$data1height <- 200
    plot.params$data2height <- 200
    
    kp <- plotKaryotype(genome = "hg19", chromosomes= paste0('chr', chrom), plot.type = 2, plot.params = plot.params)
    kpPlotRegions(kp, data = loh_seg, col="#B2182B", r0=0, r1=0.05, data.panel = 2)
    kpPlotRegions(kp, data = tai_seg, col="#2166AC", r0=0.07, r1=0.12, data.panel = 2)
    kpPlotRegions(kp, data = lst_seg, col="#9933CC", r0=0.14, r1=0.19, data.panel = 2)
    
    # plot LRR
    kpDataBackground(kp, r0 = 0, r1 = 0.55, color = "#EEFFEE", data.panel = 1)
    kpAddLabels(kp, labels="LRR", data.panel = 1, cex = 1, r0 = 0, r1 = 0.55)
    kpAxis(kp, ymin =-7, ymax = 7, r0 = 0, r1 = 0.55, numticks = 5, cex = 0.7, side = 2, data.panel = 1)
    kpPoints(kp, chr = seqz$chromosome, x = seqz$position, y = seqz$logR, col = "#AADDAA", pch=".", ymin =-7, ymax = 7, 
             r0 = 0, r1 = 0.55, data.panel = 1)
    
    # plot BAF
    kpDataBackground(kp, r0 = 0.65, r1 = 1.2, color = "#FFEEEE", data.panel = 1)
    kpAddLabels(kp, labels="BAF", data.panel = 1, cex = 1, r0 = 0.65, r1 = 1.2)
    kpAxis(kp, ymin =0, ymax = 1,r0 = 0.65, r1 = 1.2, numticks = 3, cex = 0.7, side = 2, data.panel = 1)
    kpAbline(kp, h = 0.5, col = "#DDAAAA", ymin = 0, ymax = 1, r0 = 0.65, r1 = 1.2, data.panel = 1)
    kpPoints(kp, chr = seqz$chromosome, x = seqz$position, y = seqz$Af, col = "#DDAAAA", ymin = 0, ymax = 1, data.panel = 1, r0 = 0.65, r1 = 1.2)
    kpPoints(kp, chr = seqz$chromosome, x = seqz$position, y = seqz$Bf, col = "#DDAAAA", ymin = 0, ymax = 1, data.panel = 1, r0 = 0.65, r1 = 1.2)
    
    
    legend(x = 0.04, y = 0.3,
           legend=c("LOH", "TAI", "LST"), 
           col=c("#B2182B", "#2166AC", "#9933CC"),
           pch=c(15, 15, 15))
    
    dev.off()
}
# dev.off()


