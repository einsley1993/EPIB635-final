# Date		2022-11-28
# Function	Run YL script to calculate optimal sample size, as well as futility analyses

##### Overview #####
# 0) Load packages
# 1) Running beta binomial to get posterior predictive probability of success at final analysis
# 2) Simulating the trial, no MAP
## 2.a) Writing simulation function
## 2.b) Fasle Positive Rate
## 2.c) True Positive Rate
## 2.d) Early stopping
# 3) Simulating the trial, with MAP
## 3.a) Creating function that calculates w_bar
## 3.b) False Positive Rate
## 3.c) True Positive Rate
## 3.d) Early stopping
# 4) Identifying relationship between w, TE, and FPR/TPR

##### Overview #####


# 0) Load packages
lib.loc.c <- "~/analyses/rpackages/4.0.4/cran/";
library(VGAM, lib.loc = lib.loc.c);


# 1) Running beta binomial to get posterior predictive probability of success at final analysis
# Sizes yielding 70% efficacy are nt = 30, nc = 102
pbetabinom.ab(
	117, # yf = number of cases in tx group at final analysis
	size = 264, # nf = number of cases in tx + ctl group at final
	shape1 = 1 + 30,	# a + yi, where a = 1 (noninformative distribution) and yi events in treatment group at interim = 30
	shape2 = 1 + 102, 	# b + ni - yi, where ni = total number of events at interim = 132 and yi events in treatment group at interim = 30
	log = F);
# 0.9999907
# Interpretation: At interim analysis 1, where we have 30 events in treatment group, what is the probability of observing 117 events in the treatment group at final analysis
# QUESTION: Why do we want to see 117 events at final? 117 events at final would indicate a negative study no?


# 2) Simulating the trial, no MAP
# Efficacy threshold can be specified as Un = 0.995 and U = 0.986 for final analysis
## 2.a) Writing simulation function
trial <- function(
	TE,
	ni = c(132, 264), # total number of cases at interim and final analyses
	a0 = 0, b0 = 1){
		n <- c(ni[1], diff(ni, 1));
		theta <- (1-TE)/(2-TE);
		p_eff <- NULL;
		p_success <- NULL;
		y_i <- 0;
		i <- 0;
		TE_null <- 0;
		theta_null <- (1-TE_null)/(2-TE_null);
		repeat{
			i <- i+1;

			# Simulating number of events in treatment at interim, given theta probability
			y_i <- y_i + rbinom(n = 1, size = n[i], prob = theta); 		# y_i ~ Binom(n_i, theta)

			# Updating shape parameters based on y_i
			ai <- a0 + y_i;
			bi <- b0 + ni[i] - y_i;
			
			# Efficacy: the probability of TE ≥ 0, given this posterior distribution
			peff_new <- pbeta(q = TE_null, shape1 = ai,shape2 = bi); 			
			# QUESTION: Should we use TE ≥ 0 or ≥ 0.3?
			p_eff <- c(p_eff, peff_new);

			# Success: the probability of observing 117 events in treatment at final, given events in treatment seen at this simulated interim (based on theta)
			p_success <- c(p_success, pbetabinom.ab(
						q = 117,
                                           	size = 264,
                                           	shape1 = ai,
                                           	shape2 = bi,
                                           	log = FALSE));

			# Stopping for efficacy: Is p_eff > 0.995; Stopping for futility: Is p_success < 0.05
			# We specify when to break this loop. Either stop for efficacy or stop for futility at interim, or we go on to final analysis
			# If we go onto final analysis, it may or may not have p_eff > 0.986
			if (i<=1){
				if (peff_new > 0.995 | p_success[length(p_success)] < 0.05) break;
			} else break;
		};
		out = list(p_eff, p_success);
  		return(out);
};


trial(0.2);
# [[1]]
# [1] 1
# 
# [[2]]
# [1] 0.9941065


## 2.b) False Positive Rate
fp <- NULL;
for(m in 1:10000){

	# Generate simulations under the null hypothesis (Treatment Efficacy = 0)
	p_eff <- trial(TE = 0);

	# Extracting probability of observing TE ≥ 0 (p_eff)
	p_eff_unlist <- unlist(p_eff[[1]]);

	# Indicating whether p_eff > 0.995 (threshold at interim) or 0.986 (threshold at final)
  	check <- ifelse(
		test = length(p_eff_unlist) < 2,
		yes = p_eff_unlist[length(p_eff_unlist)] > 0.995,
                no = p_eff_unlist[length(p_eff_unlist)] > 0.986);

	# We indicate whether this single trial passes threshold (True) or not (False)
	fp <- c(fp, check);
};
mean(fp);
sum(is.na(fp));
## NB: Getting 1? and only TRUE
## NB: If you get mean(fp) = 0.30. This means there is a 30% chance seeing P(TE≥0|data under null) > Threshold


## 2.c) True Positive Rate
tp <- NULL;
for(m in 1:10000){
	p_eff <- trial(TE = 0.7);
	p_eff_unlist <- unlist(p_eff[[1]]);	
	check <- ifelse(length(p_eff_unlist) < 2,
        	p_eff_unlist[length(p_eff_unlist)] > 0.995,
                p_eff_unlist[length(p_eff_unlist)] > 0.986);	
	tp <- c(tp, check);
};
mean(tp);
sum(is.na(tp));
## NB: Getting 1?


## 2.d) Early stopping
early <- NULL;
check1 <- NULL;
for(m in 1:10000){
	# Simulating a single trial under the alternative hypothesis
    	p_eff <- trial(TE = 0.7);
    	p_eff_unlist <- unlist(p_eff[[1]]);

	# Determinig whether it is shorter than 2 (it stops at the ith interim study once efficacy reached. Thus if < 2 length, means it stopped before 2nd analyis (final))
    	check <- length(p_eff_unlist) < 2;
    	early <- c(early, check);

	# Checking at which interim it stopped at specifically. Here, only one interim, so really all length(p_eff_unlist) will be one
    	check1 <- c(check1, length(p_eff_unlist)==1);
};
mean(early);
## [1] 1
mean(check1);
## [1] 1
## Here, mean(early) == mean(check1) since only one interim analysis


# 3) Simulating the trial, with MAP
## 3.a) Creating function that calculates w_bar
# posterior weight: corresponds to (1 - w_R) in Eq. (7) in Schmidli et al. (2014)
w_bar <- function(w, y, n, a, b, a0 = 1, b0 = 1){
	# QUESTION: What are these variables?
	# w = weight you initially give
	# y = events in 
	# n = total number of events
	# a = ?
	# b = ?
	f0 <- beta(a0 + y, b0 + n - y)/beta(a0, b0); 
	f <- beta(a + y, b + n - y)/beta(a, b);
	w_bar <- (1-w)*f0/((1-w)*f0 + w*f);
	return(1 - w_bar);
};

w_bar(w = 0.9, y = 20, n = 100, a = 5, b = 5);
# [1] 0.8138288


# QUESTION: w_bar(w = 0.1, y = 135, n = 132, a = 136, b = 130); w_bar may be NA. This occurs when n < y, namely n[i] < y_i, when generated in trial1


## 3.b) Writting simulation function, incorporating MAP
options(warn = -1);
trial1 <- function(
	TE, 
	w,
	ni = c(132, 264), # total no. of cases at interim and final analyses
        a0 = 1, 
	b0 = 1){
		n <- c(ni[1], diff(ni, 1));
  		theta <- (1 - TE)/(2 - TE);
 		p_eff <- NULL;
  		p_success <- NULL;
  		y_i <- 0;
		i <- 0;
		TE_null <- 0;
		theta_null <- (1-TE_null)/(2-TE_null);
		repeat{
			i <- i + 1;
			y_i <- y_i + rbinom(1, n[i], theta); # y_i ~ Binom(n_i, theta), simulating the number of cases in treatment group at interim
			ai <- a0 + y_i;
			bi <- b0 + ni[i] - y_i;

    		### Add robust MAP prior
    		# posterior weight:
		wb <- w_bar(w, y_i, n[i], ai, bi);

		# P(TE > 0|data):
		# Incorporating historical control data via yh (events at final in historical control) and nh (yh + events in historical treatment = 226 + 180)
		peff_new <- wb*pbeta(q = theta_null, shape1 = ai + 226, shape2 = bi + 226 + 180 - 226) + # shape1 = ai + yh, shape2 = bi + nh - yh
			(1-wb)*pbeta(q = theta_null, shape1 = ai, shape2 = bi);
    		p_eff = c(p_eff, peff_new);

    		p_success = c(p_success, 
			wb*pbetabinom.ab(q = 117, size = 264,
                       		shape1 = ai,
                       		shape2 = bi,
                       		log = FALSE)) +
			(1-wb)*pbetabinom.ab(q = 117, size = 264,
				shape1 = ai + 226,
				shape2 = bi + 226 + 180 - 226,
				log = FALSE);

		if (i<=1){
		# P(TE > 0|data) > 0.995;
			if (peff_new > 0.995 | p_success[length(p_success)] < 0.05) break;
			} else break ;
		};
		out = list(round(p_eff, 2), round(p_success, 3));
  		return(out);
};

trial1(TE = 0, w = 0.1);
# [[1]]
# [1] 0.93 1.00
# 
# [[2]]
# [1] 1.789 0.833





## 3.b) False Positive Rate
options(warn = -1) ;
fp <- NULL;
p_eff <- p_eff_unlist <- check <- list();
for(m in 1:10000){
  p_eff <- trial1(TE = 0, w = 0.1);
  p_eff_unlist  <- unlist(p_eff[[1]]);
  check <- ifelse(length(p_eff_unlist) < 2,
                 p_eff_unlist[length(p_eff_unlist)] > 0.995,
                 p_eff_unlist[length(p_eff_unlist)] > 0.986);
#  if(is.na(check)){
#	     break
#		 };
  fp = c(fp, check);
}; 
head(fp);
mean(fp, na.rm = T);
sum(is.na(fp));
# [1] 1644


## 3.c) True Positive Rate
options(warn = -1);
tp <- NULL
for(m in 1:10000){
	p_eff <- trial1(TE = 0.7, w = 0.5);
	p_eff_unlist <- unlist(p_eff[[1]]);
	check <- ifelse(length(p_eff_unlist) < 2,
			p_eff_unlist[length(p_eff_unlist)] > 0.995,
			p_eff_unlist[length(p_eff_unlist)] > 0.986);
	tp <- c(tp, check);
};
mean(tp, na.rm = T);
# [1] 1



## 3.d) Early stopping
# The probability that the trial stops for efficacy at the interim analysis:
options(warn = -1);
early <- NULL;
check1 <- NULL;
for(m in 1:10000){
	p_eff <- trial1(TE = 0.7, w = 0.3);
	p_eff_unlist <- unlist(p_eff[[1]]);
	check <- length(p_eff_unlist) < 2;
	early <- c(early, check);
	check1 <- c(check1, length(p_eff_unlist) == 1);
};
mean(early);
# [1] 1e-04

mean(check1);
# [1] 1e-04



# 4) Identifying relationship between w, TE, and FPR/TPR
w_change <- seq(0.1, 0.9, 0.1);


## 4.a) FPR
w_FPR <- as.data.frame(matrix(rep(0, length(w_change)*2), nrow = length(w_change), ncol = 2));
names(w_FPR) <- c("w", "FPR");
w_FPR$w <- w_change;

pb <- txtProgressBar(min = 1, max = length(w_change));
for(i in 1:length(w_change)){

	w_ix <- w_change[i];

	# Calculate fpr
	options(warn = -1) ;
	fp <- NULL;
	p_eff <- p_eff_unlist <- check <- list();
	for(m in 1:1000){
		p_eff <- trial1(TE = 0, w = w_ix);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		fp = c(fp, check);

	}; 

	# Incorporate into data.frame
	w_FPR[which(w_FPR$w == w_ix), "FPR"] <- mean(fp, na.rm = T);
	setTxtProgressBar(pb, i);
};


## 4.b) TPR
w_TPR <- as.data.frame(matrix(rep(0, length(w_change)*2), nrow = length(w_change), ncol = 2));
names(w_TPR) <- c("w", "TPR");
w_TPR$w <- w_change;

pb <- txtProgressBar(min = 1, max = length(w_change));
for(i in 1:length(w_change)){

	w_ix <- w_change[i];

	# Calculate fpr
	options(warn = -1) ;
	fp <- NULL;
	p_eff <- p_eff_unlist <- check <- list();
	for(m in 1:1000){
		p_eff <- trial1(TE = 0.7, w = w_ix);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		fp = c(fp, check);

	}; 

	# Incorporate into data.frame
	w_TPR[which(w_TPR$w == w_ix), "TPR"] <- mean(fp, na.rm = T);
	setTxtProgressBar(pb, i);
};







