'%out%' <- Negate('%in%')
subtype.colors <- data.frame(subtypes=c('AR+/NE-','ARL/NE-','AR-/NE-','AR-/NE+','AR+/NE+'),colors=RColorBrewer::brewer.pal(n=5,name = 'Dark2'))
subtype.colors$subtypes <- factor(subtype.colors$subtypes,levels=c('AR-/NE-','AR-/NE+','ARL/NE-','AR+/NE-','AR+/NE+'))
subtype.colors <- subtype.colors[order(subtype.colors$subtypes),]

match.idx = function(A, B, allow.multiple.B=F){
    if( allow.multiple.B ){
        idx.B = which(B %in% A)
        idx.A = match(B[idx.B], A)
    }
    else{
        in.both = intersect(A,B)
        idx.A = match(in.both, A)
        idx.B = match(in.both, B)
    }
    C= data.frame(idx.A, idx.B)
    if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
        stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
    C
}
alter_fun =list(
    background = function(x, y, w, h)
        grid.rect(x, y, w, h,
                  gp = gpar(fill = "white", col = "gray50", lwd = 0.4)),
    "AMP" = function(x, y, w, h) {
        grid.polygon(
            unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w),
            unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
            gp = gpar(fill = 'firebrick1', col = "black",lwd = 0.4))
    },
    DEL = function(x, y, w, h) {
        grid.polygon(
            unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w),
            unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
            gp = gpar(fill = "royalblue", col = "black",lwd = 0.4))
    },
    Frameshift = function(x, y, w, h) {
        grid.polygon(
            unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
            unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
            gp = gpar(fill = "purple", col = "black",lwd = 0.4))
    },

    "Splice donor" = function(x, y, w, h) {
        grid.polygon(
            unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
            unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
            gp = gpar(fill = "sienna3", col = "black",lwd = 0.4))
    },
    "Splice acceptor" = function(x, y, w, h) {
        grid.polygon(
            unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
            unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
            gp = gpar(fill = "orange1", col = "black",lwd = 0.4))
    },
    "Stop gained" = function(x, y, w, h) {
        grid.polygon(
            unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
            unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
            gp = gpar(fill = "seagreen2", col = "black",lwd = 0.4))
    },
    "Stop lost" = function(x, y, w, h) {
        grid.polygon(
            unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
            unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
            gp = gpar(fill = "navajowhite4", col = "black",lwd = 0.4))
    },
    "Missense" = function(x, y, w, h) {
        grid.polygon(
            unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
            unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
            gp = gpar(fill = "forestgreen", col = "black",lwd = 0.4))
    },
    biallelic = function(x, y, w, h)
        grid.points(x, y, pch = 21, size = unit(2.2, "mm"),gp=gpar(fill='black',col='white',lwd = 0.7)),
    germline = function(x, y, w, h)
        grid.points(x, y, pch = 17, size = unit(2.2, "mm"),gp=gpar(fill='white',col='hotpink2',lwd = 0.7)),
    HRD = function(x, y, w, h) {
        grid.segments(x - w*0.4, y - h*0.4, x + w*0.4, y + h*0.4, gp = gpar(col='tomato4',lwd = 1))
        grid.segments(x + w*0.4, y - h*0.4, x - w*0.4, y + h*0.4, gp = gpar(col='tomato4',lwd = 1))})
