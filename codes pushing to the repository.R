#FRANCISCO PENAGARICANO CODE FOR MULTIPLE MANHATTAN PLOT MODIFIED TO INSERT GENES INFO AT TOP WINDOWS (ssGBLUP)
#02/22/2018

# This is the new revised Manhattan Plot When the Reviewer asks us to do 
# 2019-07-23

# This is the Manhattan plot for CI 
# 2019-11-15

#### MANHATTAN PLOT FOR INTERCEPT1 ################################################
###################################################################################
rm(list = ls())
library(ggplot2)
require(ggplot2)
require(gridExtra)
require(grid)
require(cowplot)

# setting work directory
setwd("/Volumes/My Files/Final_Paper_2019_10-09/ManhattanPlot_Nov15")


# read file with snp_sol_intercept_1
SOL = read.table("snp_sol_inter_1")
colnames(SOL) <- c("Trait", "Effect", "SNP", "Chromosome", "Position", "SNPsolution", "Weight", "Variance")

## SOL  = snp_sol File output postGS filtered by SNP that passed quality control
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$Within = seq(1:nrow(SNP))


# read the windows_variance file split into chromosome_position file
# Here, we partition the chromosomes position and variances into different columns kaning it 15 
# from 9

winvar = read.table("split_windows_var_inter_1", header = F, sep = "")               
names(winvar) = c("Trait","Effect","SNPi","SNPf","tSNP","Chri","Posi","Chrf","Posf","Chr","Loc","Variance")
dim(winvar)

#number of windows to be analysed
n = 15       

# Chromosome position and location that explained largest variances in ascending order
topwin = head(winvar[order(winvar$Variance, decreasing = T),],n)  

#SNP within each window
SNP_win <- NULL      

for (i in 1:n) {
  
  chr = topwin[i,6]
  posi = topwin[i,7]
  posf = topwin[i,9]
  Window = i
  
  SNPi = cbind(SNP[SNP$Chr == chr &  SNP$Pos >= posi & SNP$Pos <= posf, ],Window)
  SNP_win = rbind(SNP_win,SNPi)
}

#order SNP by Window and Var
SNP_win = SNP_win[order(c(SNP_win$Window,SNP_win$Variance)),] 

#split SNP by Window
a = split(SNP_win, SNP_win$Window)               

#length(a[[1]][,1]) 
#number of SNPs in Window1 
#save SNP highest %var in each Window

SNP_pk = c()
for (i in 1:length(a)) {               
  SNP_pk = rbind(SNP_pk, a[[i]][1,])
}

# SNP_pk information is the most importat information . 
# We use SNP_pk information as it has Chromosome, Position, Variances , Within and Windows
g1 = "LCP1" # Chr12
g2 = "TPT1" #Chr12
g3 = "BMP5" #Chr23
g4 = "HCRTR2" #Chr23


Genes = c(g1,g2,g3,g4)
Chr = c(12,12,23,23)
Var_plot = c(1.13,1.09,1.58,1.54)      #based on y axis position
Within = c(30208,30208,50088,50088)    #within SNP_pk (Take 3rd column in the variances)
gen = data.frame(Within,Genes,Var_plot)

#SNPs must be duplicate if SNP is flagging more than one gene
SNP_gen = merge(x = SNP, y = gen, by = "Within", all = TRUE)       

#Manhattan plot
gap = 0; xx = as.numeric(); x0 = 0

for (w in 1:29) {               #29 chromosomes (Bovine)
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)              #29 chromosomes (Bovine)

color <- c(rep(c("royalblue2","olivedrab2"),14),"royalblue2")
plot1 <- ggplot(SNP, aes(x = Within, y = Var, col = factor(Chr))) + # color is based on the factor of chromosomes
  scale_color_manual(values = color) + 
  geom_point( ) + 
  geom_point(data = SNP[SNP$Within %in% Within,], aes(x = Within, y = Var), col = c("olivedrab2","royalblue2"), size = 2.0) +
  geom_text(data = SNP_gen, x = SNP_gen$Within, y = SNP_gen$Var_plot, label = SNP_gen$Genes, col = "black", fontface = "bold", size = 1.85, hjust = 0.3, vjust = -2.5) +
  ylim(0,2.0) +
  ylab("% Genetic variance")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.text.x  = element_text(size = 6)) + 
  theme(axis.text.y  = element_text(size = 8)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6))

#tiff("Intercept_1_Genes_Nov15.tiff", width = 5, height = 5, units = 'in', res = 300)
#plot_grid(plot1)
#dev.off()



##########################################################################################
######### MANHATTAN PLOT FOR SLOPE 1
#########################################################################################
#rm(list = ls())
library(ggplot2)
require(ggplot2)
require(gridExtra)
require(grid)
require(cowplot)

# setting work directory
setwd("/Volumes/My Files/Final_Paper_2019_10-09/ManhattanPlot_Nov15")


# read file with snp_sol_slope_1
SOL = read.table("snp_sol_slope_1")
colnames(SOL) <- c("Trait", "Effect", "SNP", "Chromosome", "Position", "SNPsolution", "Weight", "Variance")

## SOL  = snp_sol File output postGS filtered by SNP that passed quality control
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$Within = seq(1:nrow(SNP))


# read the windows_variance file split into chromosome_position file
# Here, we partition the chromosomes position and variances into different columns kaning it 15 
# from 9

winvar = read.table("split_windows_var_slope_1", header = F, sep = "")               
names(winvar) = c("Trait","Effect","SNPi","SNPf","tSNP","Chri","Posi","Chrf","Posf","Chr","Loc","Variance")
dim(winvar)

#number of windows to be analysed
n = 15       

# Chromosome position and location that explained largest variances in ascending order
topwin = head(winvar[order(winvar$Variance, decreasing = T),],n)  

#SNP within each window
SNP_win <- NULL      

for (i in 1:n) {
  
  chr = topwin[i,6]
  posi = topwin[i,7]
  posf = topwin[i,9]
  Window = i
  
  SNPi = cbind(SNP[SNP$Chr == chr &  SNP$Pos >= posi & SNP$Pos <= posf, ],Window)
  SNP_win = rbind(SNP_win,SNPi)
}

#order SNP by Window and Var
SNP_win = SNP_win[order(c(SNP_win$Window,SNP_win$Variance)),] 

#split SNP by Window
a = split(SNP_win, SNP_win$Window)               

#length(a[[1]][,1]) 
#number of SNPs in Window1 
#save SNP highest %var in each Window

SNP_pk = c()
for (i in 1:length(a)) {               
  SNP_pk = rbind(SNP_pk, a[[i]][1,])
}

# SNP_pk information is the most importat information . 
# We use SNP_pk information as it has Chromosome, Position, Variances , Within and Windows
g1 = "EXD2" # Chr10
g2 = "SLC10A1" #Chr10
g3 = "ADAM20" #Chr10
g4 = "HCRTR2" #Chr10
g5 = "FSD2" # Chr21
g6 = "AP3B2" # Chr21



Genes = c(g1,g2,g3,g4,g5,g6)
Chr = c(10,10,10,10,21,21)
Var_plot = c(1.09,1.06,1.03,1.00,1.02,0.99)      #based on y axis position
Within = c(26681,26681,26681,26681,47193,47193)    #within SNP_pk (Take 3rd column in the variances)
gen = data.frame(Within,Genes,Var_plot)

#SNPs must be duplicate if SNP is flagging more than one gene
SNP_gen = merge(x = SNP, y = gen, by = "Within", all = TRUE)       

#Manhattan plot
gap = 0; xx = as.numeric(); x0 = 0

for (w in 1:29) {               #29 chromosomes (Bovine)
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)              #29 chromosomes (Bovine)

color <- c(rep(c("royalblue2","olivedrab2"),14),"royalblue2")
plot2 <- ggplot(SNP, aes(x = Within, y = Var, col = factor(Chr))) + # color is based on the factor of chromosomes
  scale_color_manual(values = color) + 
  geom_point( ) + 
  geom_point(data = SNP[SNP$Within %in% Within,], aes(x = Within, y = Var), col = c("olivedrab2","blue2"), size = 2) +
  geom_text(data = SNP_gen, x = SNP_gen$Within, y = SNP_gen$Var_plot, label = SNP_gen$Genes, col = "black", fontface = "bold", size = 1.85, hjust = 0.55, vjust = -2.0) +
  ylim(0,1.5) +
  ylab("% Genetic variance")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.text.x  = element_text(size = 6)) + 
  theme(axis.text.y  = element_text(size = 8)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6))

#tiff("Slope_1_GenesJuly.tiff", width = 5, height = 5, units = 'in', res = 300)
#plot_grid(plot2)
#dev.off()

###########################################################################
# MANHATTAN PLOT FOR INTERCEPT2
###########################################################################
###
#rm(list = ls())
library(ggplot2)
require(ggplot2)
require(gridExtra)
require(grid)
require(cowplot)

# setting work directory
#setwd("/Volumes/My Files/PUBLICATION_MS_THESIS/Majorgenes_ManhattanPlot")


# read file with snp_sol_intercept_1
SOL = read.table("snp_sol_inter_2")
colnames(SOL) <- c("Trait", "Effect", "SNP", "Chromosome", "Position", "SNPsolution", "Weight", "Variance")

## SOL  = snp_sol File output postGS filtered by SNP that passed quality control
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$Within = seq(1:nrow(SNP))


# read the windows_variance file split into chromosome_position file
# Here, we partition the chromosomes position and variances into different columns kaning it 15 
# from 9

winvar = read.table("split_windows_var_inter_2", header = F, sep = "")               
names(winvar) = c("Trait","Effect","SNPi","SNPf","tSNP","Chri","Posi","Chrf","Posf","Chr","Loc","Variance")
dim(winvar)

#number of windows to be analysed
n = 15       

# Chromosome position and location that explained largest variances in ascending order
topwin = head(winvar[order(winvar$Variance, decreasing = T),],n)  

#SNP within each window
SNP_win <- NULL      

for (i in 1:n) {
  
  chr = topwin[i,6]
  posi = topwin[i,7]
  posf = topwin[i,9]
  Window = i
  
  SNPi = cbind(SNP[SNP$Chr == chr &  SNP$Pos >= posi & SNP$Pos <= posf, ],Window)
  SNP_win = rbind(SNP_win,SNPi)
}

#order SNP by Window and Var
SNP_win = SNP_win[order(c(SNP_win$Window,SNP_win$Variance)),] 

#split SNP by Window
a = split(SNP_win, SNP_win$Window)               

#length(a[[1]][,1]) 
#number of SNPs in Window1 
#save SNP highest %var in each Window

SNP_pk = c()
for (i in 1:length(a)) {               
  SNP_pk = rbind(SNP_pk, a[[i]][1,])
}

# SNP_pk information is the most importat information . 
# We use SNP_pk information as it has Chromosome, Position, Variances , Within and Windows
g1 = "LCP1" # Chr12
g2 = "TPT1" #Chr12
g3 = "BMP5" #Chr23
g4 = "HCRTR2" #Chr23


Genes = c(g1,g2,g3,g4)
Chr = c(12,12,23,23)
Var_plot = c(0.87,0.83,1.65,1.61)      #based on y axis position
Within = c(30208,30208,50087,50087)    #within SNP_pk (Take 3rd column in the variances)
gen = data.frame(Within,Genes,Var_plot)

#SNPs must be duplicate if SNP is flagging more than one gene
SNP_gen = merge(x = SNP, y = gen, by = "Within", all = TRUE)       

#Manhattan plot
gap = 0; xx = as.numeric(); x0 = 0

for (w in 1:29) {               #29 chromosomes (Bovine)
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)              #29 chromosomes (Bovine)
color <- c(rep(c("royalblue2","olivedrab2"),14),"royalblue2")
plot3 <- ggplot(SNP, aes(x = Within, y = Var, col = factor(Chr))) + # color is based on the factor of chromosomes
  scale_color_manual(values = color) + 
  geom_point( ) + 
  geom_point(data = SNP[SNP$Within %in% Within,], aes(x = Within, y = Var), col = c("olivedrab2","royalblue2"), size = 2.0) +
  geom_text(data = SNP_gen, x = SNP_gen$Within, y = SNP_gen$Var_plot, label = SNP_gen$Genes, col = "black", fontface = "bold", size = 1.85, hjust = 0.3, vjust = -2.5) +
  ylim(0,2.0) +
  ylab("% Genetic variance")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.text.x  = element_text(size = 6)) + 
  theme(axis.text.y  = element_text(size = 8)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6))

#tiff("Intercept_2_Genes_Nov15.tiff", width = 5, height = 5, units = 'in', res = 300)
#plot_grid(plot3)
#dev.off()




#tiff("Intercept_2_Genes_July.tiff", width = 5, height = 5, units = 'in', res = 300)
#plot_grid(plot1)
#dev.off()


###############################################################################
### MANHATTAN PLOT FOR SLOPE 2 #################################################
###############################################################################
#rm(list = ls())
library(ggplot2)
require(ggplot2)
require(gridExtra)
require(grid)
require(cowplot)

# setting work directory
#setwd("/Volumes/My Files/PUBLICATION_MS_THESIS/Majorgenes_ManhattanPlot")


# read file with snp_sol_slope_1
SOL = read.table("snp_sol_slope_2")
colnames(SOL) <- c("Trait", "Effect", "SNP", "Chromosome", "Position", "SNPsolution", "Weight", "Variance")

## SOL  = snp_sol File output postGS filtered by SNP that passed quality control
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$Within = seq(1:nrow(SNP))


# read the windows_variance file split into chromosome_position file
# Here, we partition the chromosomes position and variances into different columns kaning it 15 
# from 9

winvar = read.table("split_windows_var_slope_2", header = F, sep = "")               
names(winvar) = c("Trait","Effect","SNPi","SNPf","tSNP","Chri","Posi","Chrf","Posf","Chr","Loc","Variance")
dim(winvar)

#number of windows to be analysed
n = 15       

# Chromosome position and location that explained largest variances in ascending order
topwin = head(winvar[order(winvar$Variance, decreasing = T),],n)  

#SNP within each window
SNP_win <- NULL      

for (i in 1:n) {
  
  chr = topwin[i,6]
  posi = topwin[i,7]
  posf = topwin[i,9]
  Window = i
  
  SNPi = cbind(SNP[SNP$Chr == chr &  SNP$Pos >= posi & SNP$Pos <= posf, ],Window)
  SNP_win = rbind(SNP_win,SNPi)
}

#order SNP by Window and Var
SNP_win = SNP_win[order(c(SNP_win$Window,SNP_win$Variance)),] 

#split SNP by Window
a = split(SNP_win, SNP_win$Window)               

#length(a[[1]][,1]) 
#number of SNPs in Window1 
#save SNP highest %var in each Window

SNP_pk = c()
for (i in 1:length(a)) {               
  SNP_pk = rbind(SNP_pk, a[[i]][1,])
}

# SNP_pk information is the most importat information . 
# We use SNP_pk information as it has Chromosome, Position, Variances , Within and Windows
g1 = "EPAS1" # Chr11
g2 = "TAOK3" # Chr17
g3 = "NOS1" # Chr17



Genes = c(g1,g2,g3)
Chr = c(11,17,17)
Var_plot = c(1.40,0.77,0.74)      #based on y axis position
Within = c(28092,40990,40990)    #within SNP_pk (Take 3rd column in the variances)
gen = data.frame(Within,Genes,Var_plot)

#SNPs must be duplicate if SNP is flagging more than one gene
SNP_gen = merge(x = SNP, y = gen, by = "Within", all = TRUE)       

#Manhattan plot
gap = 0; xx = as.numeric(); x0 = 0

for (w in 1:29) {               #29 chromosomes (Bovine)
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)              #29 chromosomes (Bovine)
color <- c(rep(c("royalblue2","olivedrab2"),14),"royalblue2")
plot4 <- ggplot(SNP, aes(x = Within, y = Var, col = factor(Chr))) + # color is based on the factor of chromosomes
  scale_color_manual(values = color) + 
  geom_point( ) + 
  geom_point(data = SNP[SNP$Within %in% Within,], aes(x = Within, y = Var), col = c("blue2","blue2"), size = 2) +
  geom_text(data = SNP_gen, x = SNP_gen$Within, y = SNP_gen$Var_plot, label = SNP_gen$Genes, col = "black", fontface = "bold", size = 1.85, hjust = 0.55, vjust = -2.0) +
  ylim(0,1.5) +
  ylab("% Genetic variance")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.text.x  = element_text(size = 6)) + 
  theme(axis.text.y  = element_text(size = 8)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6))

#tiff("Slope_2_GenesJuly.tiff", width = 5, height = 5, units = 'in', res = 300)
#plot_grid(plot4)
#dev.off()



####################################################################################################
### MANHATTAN PLOT FOR INTERCEPT3
##################################################################################################
library(ggplot2)
require(ggplot2)
require(gridExtra)
require(grid)
require(cowplot)

# setting work directory
#setwd("/Volumes/My Files/PUBLICATION_MS_THESIS/Majorgenes_ManhattanPlot")


# read file with snp_sol_intercept_3
SOL = read.table("snp_sol_inter_3")
colnames(SOL) <- c("Trait", "Effect", "SNP", "Chromosome", "Position", "SNPsolution", "Weight", "Variance")

## SOL  = snp_sol File output postGS filtered by SNP that passed quality control
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$Within = seq(1:nrow(SNP))


# read the windows_variance file split into chromosome_position file
# Here, we partition the chromosomes position and variances into different columns kaning it 15 
# from 9

winvar = read.table("split_windows_var_inter_3", header = F, sep = "")               
names(winvar) = c("Trait","Effect","SNPi","SNPf","tSNP","Chri","Posi","Chrf","Posf","Chr","Loc","Variance")
dim(winvar)

#number of windows to be analysed
n = 15       

# Chromosome position and location that explained largest variances in ascending order
topwin = head(winvar[order(winvar$Variance, decreasing = T),],n)  

#SNP within each window
SNP_win <- NULL      

for (i in 1:n) {
  
  chr = topwin[i,6]
  posi = topwin[i,7]
  posf = topwin[i,9]
  Window = i
  
  SNPi = cbind(SNP[SNP$Chr == chr &  SNP$Pos >= posi & SNP$Pos <= posf, ],Window)
  SNP_win = rbind(SNP_win,SNPi)
}

#order SNP by Window and Var
SNP_win = SNP_win[order(c(SNP_win$Window,SNP_win$Variance)),] 

#split SNP by Window
a = split(SNP_win, SNP_win$Window)               

#length(a[[1]][,1]) 
#number of SNPs in Window1 
#save SNP highest %var in each Window

SNP_pk = c()
for (i in 1:length(a)) {               
  SNP_pk = rbind(SNP_pk, a[[i]][1,])
}

# SNP_pk information is the most importat information . 
# We use SNP_pk information as it has Chromosome, Position, Variances , Within and Windows
g1 = "LCP1" # Chr12
g2 = "TPT1" #Chr12
g3 = "BMP5" #Chr23
g4 = "HCRTR2" #Chr23

Genes = c(g1,g2,g3,g4)
Chr = c(12,12,23,23)
Var_plot = c(0.91,0.87,1.60,1.56)      #based on y axis position
Within = c(30208,30208,50087,50087)    #within SNP_pk (Take 3rd column in the variances)
gen = data.frame(Within,Genes,Var_plot)





#SNPs must be duplicate if SNP is flagging more than one gene
SNP_gen = merge(x = SNP, y = gen, by = "Within", all = TRUE)       

#Manhattan plot
gap = 0; xx = as.numeric(); x0 = 0

for (w in 1:29) {               #29 chromosomes (Bovine)
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)              #29 chromosomes (Bovine)

plot5 <- ggplot(SNP, aes(x = Within, y = Var, col = factor(Chr))) + # color is based on the factor of chromosomes
  scale_color_manual(values = color) + 
  geom_point( ) + 
  geom_point(data = SNP[SNP$Within %in% Within,], aes(x = Within, y = Var), col = c("olivedrab2","royalblue2"), size = 2.0) +
  geom_text(data = SNP_gen, x = SNP_gen$Within, y = SNP_gen$Var_plot, label = SNP_gen$Genes, col = "black", fontface = "bold", size = 1.85, hjust = 0.3, vjust = -2.5) +
  ylim(0,2.0) +
  ylab("% Genetic variance")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.text.x  = element_text(size = 6)) + 
  theme(axis.text.y  = element_text(size = 8)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6))

#tiff("Intercept_3_Genes_Nov15.tiff", width = 5, height = 5, units = 'in', res = 300)
#plot_grid(plot5)
#dev.off()
#############################################################################################################
### Manhattan Plot for slope 3
############################################################################################################
#rm(list = ls())
library(ggplot2)
require(ggplot2)
require(gridExtra)
require(grid)
require(cowplot)

# setting work directory
#setwd("/Volumes/My Files/PUBLICATION_MS_THESIS/Majorgenes_ManhattanPlot")


# read file with snp_sol_slope_3
SOL = read.table("snp_sol_slope_3")
colnames(SOL) <- c("Trait", "Effect", "SNP", "Chromosome", "Position", "SNPsolution", "Weight", "Variance")

## SOL  = snp_sol File output postGS filtered by SNP that passed quality control
SNP = data.frame(Chr = SOL$Chromosome, Pos = SOL$Position, Var = SOL$Variance)
SNP$Chr = as.numeric(as.character(SNP$Chr))
SNP$Within = seq(1:nrow(SNP))


# read the windows_variance file split into chromosome_position file
# Here, we partition the chromosomes position and variances into different columns kaning it 15 
# from 9

winvar = read.table("split_windows_var_slope_3", header = F, sep = "")               
names(winvar) = c("Trait","Effect","SNPi","SNPf","tSNP","Chri","Posi","Chrf","Posf","Chr","Loc","Variance")
dim(winvar)

#number of windows to be analysed
n = 15       

# Chromosome position and location that explained largest variances in ascending order
topwin = head(winvar[order(winvar$Variance, decreasing = T),],n)  

#SNP within each window
SNP_win <- NULL      

for (i in 1:n) {
  
  chr = topwin[i,6]
  posi = topwin[i,7]
  posf = topwin[i,9]
  Window = i
  
  SNPi = cbind(SNP[SNP$Chr == chr &  SNP$Pos >= posi & SNP$Pos <= posf, ],Window)
  SNP_win = rbind(SNP_win,SNPi)
}

#order SNP by Window and Var
SNP_win = SNP_win[order(c(SNP_win$Window,SNP_win$Variance)),] 

#split SNP by Window
a = split(SNP_win, SNP_win$Window)               

#length(a[[1]][,1]) 
#number of SNPs in Window1 
#save SNP highest %var in each Window

SNP_pk = c()
for (i in 1:length(a)) {               
  SNP_pk = rbind(SNP_pk, a[[i]][1,])
}

# SNP_pk information is the most importat information . 
# We use SNP_pk information as it has Chromosome, Position, Variances , Within and Windows
g1 = "BRWD1" # Chr1
g2 = "BMP5" #Chr23
g3 = "HCRTR2" #Chr23



Genes = c(g1,g2,g3)
Chr = c(1,23,23)
Var_plot = c(0.86,1.30,1.27)      #based on y axis position
Within = c(2969,50088,50088)    #within SNP_pk (Take 3rd column in the variances)
gen = data.frame(Within,Genes,Var_plot)

#SNPs must be duplicate if SNP is flagging more than one gene
SNP_gen = merge(x = SNP, y = gen, by = "Within", all = TRUE)       

#Manhattan plot
gap = 0; xx = as.numeric(); x0 = 0

for (w in 1:29) {               #29 chromosomes (Bovine)
  x1 = table(SNP$Chr)[[w]]
  x2 = gap + 0.5 * x1
  x3 = x0 + x2
  xx[w] = x3
  x0 = x3 + 0.5 * x1
}

name = seq(1:29)              #29 chromosomes (Bovine)

color <- c(rep(c("royalblue2","olivedrab2"),14),"royalblue2")
plot6 <- ggplot(SNP, aes(x = Within, y = Var, col = factor(Chr))) + # color is based on the factor of chromosomes
  scale_color_manual(values = color) + 
  geom_point( ) + 
  geom_point(data = SNP[SNP$Within %in% Within,], aes(x = Within, y = Var), col = c("blue2","blue2"), size = 2) +
  geom_text(data = SNP_gen, x = SNP_gen$Within, y = SNP_gen$Var_plot, label = SNP_gen$Genes, col = "black", fontface = "bold", size = 1.85, hjust = 0.55, vjust = -2.0) +
  ylim(0,1.5) +
  ylab("% Genetic variance")  + xlab("Chromosome") +
  theme(axis.title.x = element_text(size = 12)) + 
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.text.x  = element_text(size = 6)) + 
  theme(axis.text.y  = element_text(size = 8)) +
  scale_x_continuous(breaks = xx, labels = name) + theme_bw() +
  theme(legend.position = "none", panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6))

#tiff("Slope_3_GenesJuly.tiff", width = 5, height = 5, units = 'in', res = 300)
#plot_grid(plot6)
#dev.off()


##############################################################################
# how to insert 6 plots into 1 plot ##
#############################################################################
tiff("Figure_int_slo_6_14chr_20191115.tiff", width = 14, height = 12, units = 'in', res = 300)
plot_grid(plot1, plot2,plot3,plot4, plot5,plot6, align = c("hv"), nrow = 3,  
          labels = c("Parity1:General Effect", "Parity1:Thermotolerance Effect","Parity2:General Effect","Parity2:Thermotolerance Effect","Parity3:General Effect","Parity3:Thermotolerance Effect"), label_size= 10,hjust = -0.50, vjust=2.5, label_colour = "darkgreen")
dev.off()