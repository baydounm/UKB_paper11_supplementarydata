#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("readstata13")

library(ggplot2)
library(ggrepel)

# setwd("D:\\DATA")
# fnm = 'VOLCANOPLOT_DATASET.dta'

fnmCSV = '/zdsk/Manuscripts/Baydoun/UK-BioBank/2023-12-20—11-ProteomicBiomarkers/Volcano/LE8_PROTEOME_VOLCANOPLOTDATA.csv'

datObj = read.csv(fnmCSV, as.is=T)
nrow(datObj)
summary(datObj)
head(datObj)

pltDat = within(datObj, {
	lblOutcome	= sub('^ z', '', zOutcome, perl=T)
# 	p		= pmax(p, p, 10^-600)
	})

summary(pltDat)
head(pltDat)

options(ggrepel.max.overlaps = Inf)	# avoids error message about colliding labels

Volcano = function(volDat) {
	ggObj =	ggplot(data=volDat, aes(x=estimate, y=mlog10p, label=lblOutcome)) +
		geom_point(shape=20, size=1) +
                geom_point(shape=20, size=3, col='red', data=subset(volDat, selected_final==1)) +
		geom_point(shape=20, size=3, col='orange', data=subset(volDat, selected_final==1 & signif==1 & estimate>0.20)) +
                geom_point(shape=20, size=3, col='blue', data=subset(volDat, selected_final==1 & signif==1 & estimate< -0.20)) +
# 		xlim(-.45, .45) + ylim(0, 250) +
		xlim(-.45, .45) + ylim(0, 540) +
		geom_text_repel(size=3, col='red', data=subset(volDat, selected_final==1 & estimate>0.20), aes(label=lblOutcome)) +
                geom_text_repel(size=3, col='blue', data=subset(volDat, selected_final==1 & estimate< -0.20), aes(label=lblOutcome))
                labs(x='Effect size') +
		theme_minimal() +
		geom_vline(xintercept=-.20, col="gray", linetype=2) +
		geom_vline(xintercept=.20,  col="gray", linetype=2) +
		geom_hline(yintercept=-log10(0.05/1463), col="gray", linetype=2)

	ggObj
	}

(pltObj = Volcano(pltDat))

dirName = '//zdsk/Manuscripts/Baydoun/UK-BioBank/2023-12-20—11-ProteomicBiomarkers/Volcano/'
imgName = '2023-12-21-UKB11-Volcano'
pdfFile = paste0(dirName, imgName, '.pdf')
pngFile = paste0(dirName, imgName, '.png')

pdf(pdfFile)
pltObj
jnk = dev.off()

png(pngFile)
pltObj
jnk = dev.off()
