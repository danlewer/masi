# This script calculates the MASI and YLLI indicators for England, 2002-2016, using publicly available data.
# It is written for R version 3.5.1, using the 'data.table' package. You may need to use the command "install.packages('data.table')".
# It shows the general method behind the results in the article:

### Lewer D, Jayatunga W, Aldridge RW, Edge C, Marmot M, Story A, Hayward A. Premature mortality
### attributable to socioeconomic inequality in England between 2003 and 2018: an observational 
### study. 2019. Lancet Public Health. DOI:10.1016/ S2468-2667(19)30219-1

# The mortality and population data is available at: 
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/adhocs/007710numberofdeathsandpopulationsindeprivationdecileareasbysexandsingleyearofageenglandandwalesregisteredyears2001to2016
# For ease, we have reformatted the data and uploaded it to a public repository, where it can be read directly using the code provided.
# Results of this code vary slightly from results published in the article due to small differences between pubicly available data and microdata used in the article, and differences in the time period studied.

library(data.table)

#-----------------------------------
# Read population and mortality data
#-----------------------------------

pop <- read.csv('https://raw.githubusercontent.com/danlewer/masi/main/pop.csv')
mort <- read.csv('https://raw.githubusercontent.com/danlewer/masi/main/mort.csv')
setDT(pop); setDT(mort)

#-------------------------------------------------------
# Combine datasets and group into age groups and periods
#-------------------------------------------------------

# combine

pop <- melt(pop, measure.vars = paste0('X', as.character(0:74)), variable.name = 'age', value.name = 'pop', variable.factor = F)
mort <- melt(mort, measure.vars =  paste0('X', as.character(0:74)), variable.name = 'age', value.name = 'deaths', variable.factor = F)
d <- mort[pop, on = c('year', 'imd10', 'sex', 'age')]

# age groups

d[, age := gsub('X', '', age)]
d[, age := as.numeric(age)]
d[, age_group := findInterval(age, c(0, 1, seq(5, 70, 5)))]
ages <- c('0', '1-4', paste0(seq(5, 70, 5), '-', seq(9, 74, 5)))
d[, age_group := factor(age_group, 1:16, ages)]

# period

d <- d[year != 2001]
d[, period := findInterval(year, c(2002, 2005, 2008, 2011, 2014))]

#---------------
# Calculate MATI
#---------------

# aggregated dataset

d_mati <- d[, .(deaths = sum(deaths), pop = sum(pop)), c('imd10', 'sex', 'age_group', 'period')]

# reference rate

ref <- d_mati[imd10 == 10]
ref[, ref_rate := deaths / pop * 100000]
ref <- ref[,c('sex', 'age_group', 'period', 'ref_rate'), with = F]

# calculated expected and attributable

d_mati <- ref[d_mati, on = c('sex', 'age_group', 'period')]
d_mati[, expected := pop * ref_rate / 100000]
d_mati[, attributable := deaths - expected]

# MASI

sum(d_mati$deaths) # 2.4m deaths
sum(d_mati$attributable) # 859k attributable
sum(d_mati$attributable) / sum(d_mati$deaths) # MATI = 36.3%

# MATI by IMD

mati_by_imd <- aggregate(cbind(expected, attributable) ~ imd10, d_mati, sum)

#-----------------------------------
# Create life tables and report YLLI
#-----------------------------------

# function to calculate total YLLs using a life table

total_yll <- function(mx_, cohort = 100000) { # mx_ is vector of mortality rates by age
  mx <- c(mx_, 1)
  n <- length(mx)
  qx <- 2 * mx / (2 + mx)
  lx <- c(1, cumprod(1 - qx))[seq_len(n)] * cohort
  dx <- -c(diff(lx), lx[n] * qx[n])
  dx <- dx[1:(n-1)]
  yll <- dx * (((n-1):1) - 0.5) # assume 74.5 years lost for death age 0; 73.5 years age 1, etc.
  return(sum(yll))
}

# aggregated dataset

d_yli <- aggregate(cbind(deaths, pop) ~ imd10 + sex + age, d, sum)

# mortality rates by age, sex and IMD10

d_yli$mx <- d_yli$deaths / d_yli$pop
d_yli <- split(d_yli$mx, f = list(d_yli$imd10, d_yli$sex))

# YLLs by age, sex and IMD10

d_yli <- sapply(d_yli, total_yll)

# add male and female together, and calculate attributable YLLs, and express per person

d_yli <- d_yli[1:10] + d_yli[11:20]
names(d_yli) <- 1:10
d_yli <- rbind(total = d_yli, expected = d_yli[10])
d_yli <- rbind(d_yli, attributable = d_yli[1,] - d_yli[2,]) / 200000

sum(d_yli[3,]) / sum(d_yli[1,]) # 37% of YLLs attributable to inequality
sum(d_yli[3,]) / 10 # 1.28 YLLs attributable to inequality per person

#---------
# Barplots
#---------

par(mfrow = c(1, 2))

barplot(t(as.matrix(mati_by_imd[10:1,-1])), col = c("#FB9A99", "#E31A1C"), space = 0, ylab = 'Premature deaths, England, 2002-2016', xlab = 'IMD', main = 'MATI')
barplot(d_yli[2:3,10:1], col = c("#FB9A99", "#E31A1C"), space = 0, ylab = 'Mean YLLs per person', xlab = 'IMD', main = 'YLLI')

#------------------------------------------
# Monte-Carlo confidence intervals for MATI
#------------------------------------------

# General function

mc_paf <- function(dat, mod_var, ref, adj_vars, time = 'pop', event = 'deaths', N = 10000, point.estimate = mean, level = 0.95) {
  
  # deaths and mortality rates
  rd <- sapply(dat[,event], function(x) rpois(N, x))
  rd_rate <- t(t(rd) / dat[,time])
  
  # which are the relevant reference groups?
  refn <- dat[,c(mod_var, adj_vars)]
  refn$nc1 <- seq_len(nrow(refn))
  refn_ref <- refn[refn[,mod_var] == ref,]
  refn_ref[,mod_var] <- NULL
  names(refn_ref) <- c(adj_vars, 'nc2')
  refn <- merge(refn, refn_ref, by = adj_vars)
  refn <- refn$nc2[order(refn$nc1)]
  
  # expected
  expected <- t(t(rd_rate[,refn]) * dat[,time])
  
  # attributable
  attributable <- rd - expected
  attributable <- attributable[,dat[,mod_var] != ref] # remove those in reference category
  attributable <- rowSums(attributable)
  
  # PAF and results
  pafs <- attributable / rowSums(rd)
  point <- point.estimate(pafs)
  c(point, quantile(pafs, c((1-level)/2, 1-(1-level)/2)))
  
}

# Estimate 95% confidence intervals

d_mati$imd2 <- ifelse(d_mati$imd10 == 10, 'high', 'low') # just for efficiency. Can also use imd10 with '10' as the reference group
d_mati2 <- aggregate(cbind(deaths, pop) ~ sex + age_group + period + imd2, d_mati, sum)

mc_paf(d_mati2, mod_var = 'imd2', ref = 'high', adj_vars = c('sex', 'age_group', 'period'))
# MATI = 36.3 (95% CI 36.0-36.6)
