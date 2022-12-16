# Date		2022-11-27
# Function	Grid search to identify events in treatment and control totalling efficacy >70%
# Author	Joan Miguel Romero




###### Overview ######
# 0) Load packages
# 1) Function to calculate efficacy based on cases in treatment and cases in control
# 2) Grid search to idenitfy which permutations yield 70% efficacy
# 3) Identifying which splits give you efficacy cutoff
# 4) Outputting results
###### Overview ######

# Summary
# Anticipate: 	117 events in treatment, 147 events in control, at final analysis (n = 2000)
# Reasoning:	117 = (180/3071)*2000, where 180 were number of events and total number of patients, respectively, in treatment group at final analysis
#		147 = (226/3076)*2000, where 226 were number of events and total number of patients, respectively, in control group at final analysis




# 0) Load packages
lib.loc.c <- "~/analyses/rpackages/4.0.4/cran/";
library(gridExtra, lib.loc = lib.loc.c);                 


# 1) Function to calculate efficacy based on cases in treatment and cases in control
efficacy <- function(nt, nc){
	n <- nt+nc;
	p.t <- nt/n;
	p.c <- nc/n;
	theta <- p.t/(p.t+p.c);
	eff <- (2*theta-1)/(theta-1);
	return(eff);
};


# 2) Grid search to idenitfy which permutations yield 70% efficacy
nt.nums <- seq(1, 99, 1);
ratios <- data.frame(nt = as.numeric(nt.nums), nc = as.numeric(99-nt.nums));
eff.results <- rep(0, nrow(ratios));

ratios.eff <- cbind(ratios, eff.results);
for(i in 1:nrow(ratios)){
	ratios.eff[i, 3] <- efficacy(ratios.eff[i, 1], ratios.eff[i, 2]);
};


# 3) Identifying which splits give you efficacy cutoff
eff.70 <- ratios.eff[which(ratios.eff$eff.results > 0.7),]
eff.70[order(eff.70$eff.results)[1], ];

# 4) Outputing results
pdf(file = "~/courses/epib635/graphs/project_4_efficacy.pdf", height = 40, width = 7);
grid.table(ratios.eff);
dev.off();


