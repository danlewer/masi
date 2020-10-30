library(data.table)
library(stringi)
library(ggthemes)

#  ------------------------------
#  function to draw a pie section
#  ==============================

circular_segment <- function(radius, deg1 = 0, deg2 = 180, centre = c(0, 0), n = 1000, COL = 'red', BORDER = NA, LWD = 0.75, type = 'segment') {
  angles <- seq(deg1, deg2, length.out = n) / (180 / pi)
  x <- sin(angles) * radius + centre[1]
  y <- cos(angles) * radius + centre[2]
  if (type == 'segment') {
    polygon(x, y, col = COL, border = BORDER, lwd = LWD)
  } else {
    polygon(c(centre[1], x), c(centre[2], y), col = COL, border = BORDER, lwd = LWD)
  }
}

#  --------------------
#  EXPLODOGRAM function
#  ====================

explodogram2 <- function(d1, d2, 
                         c1 = c(0, 0), # co-ordinates of centre circle
                         start_angle = 200, # starting angle for level 1
                         covers1 = 300, # degrees over which level 1 is arranged
                         ANGLES = NA,
                         disp = 5, # length of lines between centre circle and level 1. single number, or vector of numbers giving different lengths
                         disp2 = 2, # length of lines for level 2
                         main_rad = 2, # radius of centre circle. all other circles calculated from this
                         col_att = "#BDD5EA", # colour of attributable. Can be vector
                         col_exp = "#577399", # colour of expected. Can be vector
                         col_all = c("#f47942", "#fbb04e"), # colour of centre pie
                         col_line = "#fbb04e", # line between centre of pie and 1st level circles
                         col_line2 = "#fbb04e", # lines between 1st level circles and 2nd level circles
                         XLIM = c(-10, 10), # arbitrary horizontal size of plot (see 'c1' and 'main_rad')
                         YLIM = c(-10, 10), # arbitrary vertical size of plot
                         l1label = 'chapter', # name of top level variable in your data
                         l2label = 'major', # name of 2nd level variable in your data
                         labPos = 3, # position of labels
                         CEX = 1, # size of labels
                         CENTRE_TITLE = 'ALL', # text in the middle of the pie
                         MAIN = '', # title above the chart
                         CEX.MAIN = 2, # size of text in the middle of the pie
                         l1alignx, # adjustments to label alignments
                         l1aligny,
                         l2alignx,
                         l2aligny
                         
) {
  
  nc1 <- nrow(d1)
  col_att <- if (length(col_att) == 1L) rep(col_att, nc1) else col_att
  col_exp <- if (length(col_exp) == 1L) rep(col_exp, nc1) else col_exp
  col_line <- if (length(col_line) == 1L) rep(col_line, nc1) else col_line
  col_line2 <- if (length(col_line2) == 1L) rep(col_line2, nc1) else col_line2
  
  # location of centre of exploded circle given angle and distance from centre
  
  f <- function(angles, hyp = disp, start = c(0, 0), out = 'LIST') { # angles in degrees
    a <- abs((angles %% 180) - 90)
    a <- a * (pi / 180)
    y <- sin(a) * hyp
    x <- cos(a) * hyp
    b <- angles %% 360
    x <- ifelse(b > 180, -x, x)
    y <- ifelse(b > 90 & b < 270, -y, y) 
    x <- x + start[1]
    y <- y + start[2] 
    if (out == 'LIST') return(list(x, y)) else return(c(x, y))
  }
  
  # calculate angles and centre points: level 1
  
  if(is.na(ANGLES[1])) {
    angles <- seq(start_angle, start_angle + covers1, length.out = nc1)
  } else {
    angles <- ANGLES
  }
  sc <- mapply(f, angles, hyp = disp, out = 'CONC')
  
  # calculate angles and centre points: level 2
  
  ang2_even <- c(-20, 20, -60, 60, -100, 100, -140, 140) # offset degrees for level 2
  ang2_odd <- c(0, -40, 40, -80, 80, -120, 120, -160, 160)
  
  nl2 <- sapply(d2, nrow) # number of level-2 explosions
  f2 <- function(angle, nl) {
    ang2 <- if((nl %% 2) == 0) ang2_even else ang2_odd
    (angle %% 360) + ang2[seq_len(nl)]
  }
  angles2 <- mapply(f2, angle = angles, nl = nl2)
  
  sc1 <- split(t(sc), seq_len(nc1))
  sc2 <- mapply(f, angles = angles2, hyp = disp2, start = sc1, out = 'LIST', SIMPLIFY = F)
  
  # circle radii: level 1
  
  mrii <- colSums(d1[,c('observed', 'expected')])
  r1 <- sqrt(mrii[1] / pi)
  sf <- main_rad / r1
  r1 <- r1 * sf
  rii <- sqrt(d1[,c('observed', 'expected')] / pi) * sf
  
  # circle radii: level 2
  
  l2obs <- lapply(d2, function(x) sqrt(x$observed / pi) * sf)
  l2exp <- lapply(d2, function(x) sqrt(x$expected / pi) * sf)
  
  # base plot
  
  plot(0, 0, xlim = XLIM, ylim = YLIM, type = 'n', axes = F, xlab = NA, ylab = NA, asp = 1)
  title(main = MAIN, cex.main = CEX * CEX.MAIN)
  
  # draw lines
  
  segments(0, 0, sc[1,], sc[2,], col = col_line, lwd = 8)
  
  for(i in seq_len(nc1)[nl2 > 1]) {
    mapply(segments, x0 = sc[1,i], y0 = sc[2,i], x1 = sc2[[i]][1], y1 = sc2[[i]][2], col = col_line2[i], lwd = 3)
  }
  
  # draw circles: centre
  
  symbols(x = c1[1], y = c1[2], circles = r1, add = T, inches = F, bg = col_all[2], fg = 'white')
  seg_ang0 <- (mrii[1] - mrii[2]) / mrii[1] * 360
  circular_segment(r1, 0, seg_ang0, COL = col_all[1], type = 'pie', BORDER = 'white')
  
  # draw circles: level 1
  
  mati1 <- d1[,c('observed', 'expected')]
  mati1 <- (mati1$observed - mati1$expected) / mati1$observed
  mati1 <- ifelse(mati1 < 0, 0, mati1)
  
  symbols(x = sc[1,], y = sc[2,], circles = rii[,1], add = T, inches = F, bg = col_att, fg = 'white')
  seg_ang1 <- mati1 * 360
  for(i in seq_len(nc1)) {
    circular_segment(rii[i,1], 0, seg_ang1[i], centre = sc1[[i]], COL = col_exp[i], type = 'pie', BORDER = 'white')
  }
  
  # draw circles: level 2
  
  fsl <- function(radius, x, y, prop, COL = 'red') {
    ang <- prop * 360
    circular_segment(radius, 0, ang, c(x, y), COL = COL, type = 'pie', BORDER = 'white')
  }
  
  mati2 <- lapply(d2, function(x) (x$observed - x$expected) / x$observed)
  mati2 <- lapply(mati2, function(x) ifelse(x < 0, 0, x))
  
  for(i in seq_len(nc1)[nl2 > 1]) {
    symbols(sc2[[i]][[1]], sc2[[i]][[2]], circles = l2obs[[i]], add = T, inches = F, bg = col_att[i], fg = 'white')
    mapply(fsl, radius = l2obs[[i]], x = sc2[[i]][[1]], y = sc2[[i]][[2]], prop = mati2[[i]], COL = col_exp[i])
  }
  
  # labels
  
  text(c1[1], c1[2], CENTRE_TITLE, font = 2, cex = CEX * 1.5)
  
  text(sc[1,] + l1alignx, sc[2,] + l1aligny, d1[,l1label], cex = CEX, pos = labPos, font = 2)
  
  sc2 <- lapply(sc2, function(x) do.call('rbind', x))
  # sc2 <- lapply(sc2, function(x) if (ncol(x) > 8) x[,1:8] else x)
  
  for(i in seq_len(nc1)[nl2 > 1]) {
    text(sc2[[i]][1,] + l2alignx[[i]], sc2[[i]][2,] + l2aligny[[i]], d2[[i]][,l2label], cex = CEX, pos = labPos)
  }
  
}

#  -----------
#  format data
#  ===========

d <- read.csv('https://raw.githubusercontent.com/danlewer/masi/main/pie_chart_data.csv')
setDT(d)
d[, chapter2 := stri_replace_all_fixed(chapter2, ' ', '\n')]
d[, medium3 := stri_replace_all_fixed(medium3, ' ', '\n')]
d[, chapter2 := stri_replace_all_fixed(chapter2, '\n&', ' &')]
d[, medium3 := stri_replace_all_fixed(medium3, '\n&', ' &')]
d[, medium3 := stri_replace_all_fixed(medium3, '\n[', ' [')]
d$medium3[d$medium3 == 'Other\nheart\ndisease'] <- 'Other heart\ndisease'
d$medium3[d$medium3 == 'Veins &\nlymph\nnodes'] <- 'Veins &\nlymph nodes'
d$medium3[d$medium3 == 'Pulmonary\nheart\ndisease'] <- 'Pulmonary\nheart disease'
d$medium3[d$medium3 == 'Psychoactive\nsubstance\nuse'] <- 'Alcohol\n& drugs'
d$medium3[d$medium3 == 'Events\nof\nundetermined\nintent'] <- 'Events of\nundetermined\nintent'
d$medium3[d$medium3 == 'Stomach [1]'] <- 'Stomach\n[1]'
d1 <- aggregate(cbind(observed, expected) ~ chapter2, d, sum)
d1 <- d1[order(d1$observed, decreasing = T),]
d2 <- aggregate(cbind(observed, expected) ~ chapter2 + medium3, d, sum)
d2 <- d2[order(d2$observed, decreasing = T),]
d2 <- split(d2, f = d2$chapter)
d2 <- d2[d1$chapter]
d2 <- lapply(d2, function(x) x[x$observed > 10000,])
d2 <- lapply(d2, function(x) if(nrow(x) > 8L) x[1:8,] else x)

#  ------------
#  create chart
#  ============

# colours
cols <- tableau_color_pal('Miller Stone')(11)
cols <- cols[c(4, 3, 6, 7)]

#postscript("plot.eps",
#           horizontal = F, onefile = T, paper = "special", height = 12, width = 12, family = 'serif')

par(lheight = 0.8, xpd = T, oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))

explodogram2(d1, d2,
             #ANGLES = c(200, 245, 280, 310, 340, 370, 400, 430, 460, 485, 510),
             ANGLES = c(200, 245, 283, 305, 340, 370, 400, 430, 460, 485, 510),
             disp = c(5.5, 9, 5, 8, 4.5, 7.5, 7, 5, 5, 4, 4), 
             disp2 = c(2.7, 2.5, 2, 1.7, 1.7, 1.7, 2, 1, 1.3, 1, 1),
             main_rad = 3,
             XLIM = c(-14, 9),
             YLIM = c(-8, 8),
             labPos = NULL,
             l1label = 'chapter2',
             l2label = 'medium3',
             col_exp = cols[1],
             col_att = cols[2],
             col_all = cols[3:4],
             col_line = cols[4],
             col_line2 = cols[4],
             CEX = 1.1,
             CEX.MAIN = 0.8,
             CENTRE_TITLE = 'Premature\nmortality\n(2,456,285)',
             l1alignx = c(0, 0, 0, 0, 0, 0, 0, -1.5, 0, 1.1, 1.2),
             l1aligny = c(0, 0, 0, 0, 0, 0.9, -1, 0, 0.8, 0, 0),
             l2alignx = list(
               c(0, 0, 0, 0, 1.7, 0, 0, 0),
               c(-1.9, 0, -1.4, 0, -1.3, 0, 1.4),
               c(-1.3, -1.3, -0.9),
               c(-1.4, -1.2, 0),
               c(0, -1, 1, 0, 0),
               NA,
               c(1, 0.8, 0, 1.4, 0, 1.6),
               c(1, 0.7),
               c(1, 1),
               NA,
               NA),
             l2aligny = list(
               c(-1.3, -1.3, -0.9, -0.8, 0, 1.1, 0.6, 1),
               c(0, -0.9, 0, -0.6, 0, 0.7, 0),
               c(0, 0, 0),
               c(0, 0, 0.9),
               c(1.1, 0, 0, -0.6, -0.4),
               NA,
               c(0, 0, 0.8, 0, -0.5, 0),
               c(0, 0),
               c(0, 0),
               NA,
               NA)
)

# dev.off()
