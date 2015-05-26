#### plot2up.r, version 0.0.0a, 23 August 2013

args <- commandArgs(trailingOnly = TRUE)
plot2up.prefix = args[1]
# value of SCALE parameter used as pixel argument for high-resolution pngs:
s_factor = as.numeric(args[2])
numsites = as.numeric(args[3])
### in retrospect passing the numsites argument is superfluous 
### as it could be read from png file dimensions

library(ape)
library(png)

plot_tree = function(plot2up.tree, plot2up.ylim) {

    print(paste("tree has", length(plot2up.tree$tip.label), "taxa"), quote=F)

    treeout = plot(plot2up.tree, font=1, no.margin=F, y.lim=plot2up.ylim)
    plot2up.xlim=treeout$x.lim
    add.scale.bar(plot2up.xlim[2]*0.85, 1)

    plot2up.xlim[2]
}

plot_png = function(pngfile, plot2up.ylim, scale_factor, is_aa) {

    plot2up.pngraster = readPNG(pngfile, native=T, info=T)
    currdim<-dim(plot2up.pngraster)

    aa_scale = ifelse(is_aa, 3, 1)
    plot2up.xlim = c(1/2, -1/2+aa_scale*(currdim[2]/scale_factor))

    print(paste("pixel dimensions are w", currdim[2]/scale_factor, 
	    "x h", currdim[1]/scale_factor), quote=F)

    print(paste("scaled to w", (plot2up.xlim[2]+1/2)-(plot2up.xlim[1]-1/2), 
	    "x h", (plot2up.ylim[2]+1/2)-(plot2up.ylim[1]-1/2)), quot=F)

    rasterImage(plot2up.pngraster, plot2up.xlim[1]-1/2, plot2up.ylim[1]-1/2, plot2up.xlim[2]+1/2, plot2up.ylim[2]+1/2, interpolate=F)
}

# Start here
outfile = paste0(plot2up.prefix, '.pdf')
pdf(outfile, paper='USr', width=11, height=8.5)

layout(matrix(c(1:2), ncol=2, byrow=T))

par(ljoin=1, lend=1)
par(cex.axis=4/5, cex.lab=4/5)
par(xaxt='n', yaxt='n', tcl=-1/5)
par(oma=c(0, 0, 0, 0))
par(mar=c(0, 0, 0, 0), adj=0)
par(mgp=c(0, 0, 0))
par(xaxs='r', yaxs='r')

treefile = paste0(plot2up.prefix, '.rtree')
plot2up.tree = read.tree(treefile)
plot2up.xlim=c(1,numsites)
plot2up.ylim=c(1,length(plot2up.tree$tip.label))

# Invoke plot to establish coordinates but show nothing:
plot(0, 0, type='n', frame.plot=F, xlab='', ylab='', 
    xlim=plot2up.xlim, ylim=plot2up.ylim)

plot_png(paste0(plot2up.prefix, ".png"), plot2up.ylim, s_factor, T)

# Uncomment to add amino-acid png:
#plot_png(paste0(plot2up.prefix, "-aa.png"), plot2up.ylim, s_factor, T)

par(xaxt='s') # Show the x-axis
axis(1, pos=1/2)

# To avoid using ape's (broken) ladderize function, invert y-axis coordinates.
# Otherwise the tree would be upside-down.
xmax = plot_tree(plot2up.tree, rev(plot2up.ylim))
# xmax gives the width of the tree as rendered (min. value is 0)

print(paste("tree width =", xmax), quote=F)

# Add title to plot:
#mtext(plot2up.prefix, 3, line=-1/2, outer=T)

dev.off()
embedFonts(outfile)
system(paste("open", outfile))
