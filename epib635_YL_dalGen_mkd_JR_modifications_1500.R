# Date		2022-12-12
# Function	Run YL script to calculate optimal sample size, as well as futility analyses, using values represeting 1500 total cases

##### Overview #####
# 0) Load packages
# 1) Running beta binomial to get posterior predictive probability of success at final analysis (YL)
# 2) Simulating the trial, no MAP
## 2.a) Writing simulation function (YL)
## 2.b) Fasle Positive Rate ((YL)
## 2.c) True Positive Rate (YL)
## 2.d) Early stopping (YL)
# 3) Simulating the trial, with MAP
## 3.a) Creating function that calculates w_bar (YL)
## 3.b) False Positive Rate (YL)
## 3.c) True Positive Rate (YL)
## 3.d) Early stopping (YL)
# 4) Identifying relationship between w, TE, and FPR/TPR
## 4.a) Figure 1: FPR vs w, setting TE = 0.0 (JR)
## 4.b) Figure 2: TPR vs w, setting TE = 0.7 (JR)
## 4.c) Figure 3: TPR vs TE, setting w = 0.7 (JR)
## 4.d) Plotting (JR)
# 5) Borrowing for both control and treatment arms
## 5.a) Updating w_bar function to include treatment data (YL)
## 5.b) False Positive Rate (YL)
## 5.c) True Positive Rate (YL)
# 6) Identifying relationship between w, TE, and FPR/TPR using historical data from both arms
## 6.a) Figure 4: FPR vs w, setting TE = 0.0 (JR)
## 6.b) Figure 5: TPR vs w, setting TE = 0.7 (JR)
## 6.c) Figure 6: TPR vs TE, setting w = 0.7 (JR)
## 6.d) Figure 7: TPR vs TE, setting w = 0.5 (JR)
## 6.d) Plotting (JR)
# 7) Save image

##### Overview #####


# 0) Load packages
lib.loc.c <- "~/analyses/rpackages/4.0.4/cran/";
library(VGAM, lib.loc = lib.loc.c);
library(withr, lib.loc = lib.loc.c);		# Dependency of ggplot2
library(labeling, lib.loc = lib.loc.c);		# Dependency of ggplot2
library(farver, lib.loc = lib.loc.c);		# Dependency of ggplot2
library(digest, lib.loc = lib.loc.c);		# Dependency of ggplot2
library(ggplot2, lib.loc = lib.loc.c);


# 1) Running beta binomial to get posterior predictive probability of success at final analysis (YL)
# Sizes yielding 70% efficacy are nt = 22, nc = 77
pbetabinom.ab(
	88, # yf = number of cases in tx group at final analysis
	size = 198, # nf = number of cases in tx + ctl group at final
	shape1 = 1 + 22,	# a + yi, where a = 1 (noninformative distribution) and yi events in treatment group at interim = 22
	shape2 = 1 + 77, 	# b + ni - yi, where ni = total number of events at interim = 99 and yi events in treatment group at interim = 22
	log = F);
# 0.9999907
# Interpretation: At interim analysis 1, where we have 22 events in treatment group, what is the probability of observing 88 events in the treatment group at final analysis


# 2) Simulating the trial, no MAP (YL)
# Efficacy threshold can be specified as Un = 0.995 and U = 0.986 for final analysis
## 2.a) Writing simulation function
trial <- function(
	TE,
	ni = c(99, 198), # total number of cases at interim and final analyses
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
			p_eff <- c(p_eff, peff_new);

			# Success: the probability of observing 88 events in treatment at final, given events in treatment seen at this simulated interim (based on theta)
			p_success <- c(p_success, pbetabinom.ab(
						q = 88,
                                           	size = 198,
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


#trial(0.2);
# [[1]]
# [1] 1
# 
# [[2]]
# [1] 0.9941065


## 2.b) False Positive Rate (YL)
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
#mean(fp);
#sum(is.na(fp));
## NB: If you get mean(fp) = 0.30. This means there is a 30% chance seeing P(TE≥0|data under null) > Threshold


## 2.c) True Positive Rate (YL)
tp <- NULL;
for(m in 1:10000){
	p_eff <- trial(TE = 0.7);
	p_eff_unlist <- unlist(p_eff[[1]]);	
	check <- ifelse(length(p_eff_unlist) < 2,
        	p_eff_unlist[length(p_eff_unlist)] > 0.995,
                p_eff_unlist[length(p_eff_unlist)] > 0.986);	
	tp <- c(tp, check);
};
#mean(tp);
#sum(is.na(tp));


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
#mean(early);
## [1] 1
#mean(check1);
## [1] 1
## Here, mean(early) == mean(check1) since only one interim analysis


# 3) Simulating the trial, with MAP
## 3.a) Creating function that calculates w_bar (YL)
# posterior weight: corresponds to (1 - w_R) in Eq. (7) in Schmidli et al. (2014)
w_bar <- function(w, y, n, a, b, a0 = 1, b0 = 1){
	f0 <- beta(a0 + y, b0 + n - y)/beta(a0, b0); 
	f <- beta(a + y, b + n - y)/beta(a, b);
	w_bar <- (1-w)*f0/((1-w)*f0 + w*f);
	return(1 - w_bar);
};

w_bar(w = 0.9, y = 20, n = 100, a = 5, b = 5);
# [1] 0.8138288



## 3.b) Writting simulation function, incorporating MAP (YL)
options(warn = -1);
trial1 <- function(
	TE, 
	w,
	ni = c(99, 198), # total no. of cases at interim and final analyses
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
			wb*pbetabinom.ab(q = 88, size = 198,
                       		shape1 = ai,
                       		shape2 = bi,
                       		log = FALSE)) +
			(1-wb)*pbetabinom.ab(q = 88, size = 198,
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

#trial1(TE = 0, w = 0.1);
# [[1]]
# [1] 0.93 1.00
# 
# [[2]]
# [1] 1.789 0.833





## 3.b) False Positive Rate (YL)
options(warn = -1) ;
fp <- NULL;
for(m in 1:100){
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
#head(fp);
#mean(fp, na.rm = T);
#sum(is.na(fp));
# [1] 1644


## 3.c) True Positive Rate (YL)
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
#mean(tp, na.rm = T);
# [1] 1



## 3.d) Early stopping (YL)
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
#mean(early);
# [1] 1e-04

#mean(check1);
# [1] 1e-04


start.time <- Sys.time();
# 4) Identifying relationship between w, TE, and FPR/TPR
w_change <- seq(0, 1, 0.01);
te_change <- seq(0.3, 0.9, 0.01);


## 4.a) Figure 1: FPR vs w, setting TE = 0.0 (JR)
w_FPR <- as.data.frame(matrix(rep(0, length(w_change)*2), nrow = length(w_change), ncol = 2));
names(w_FPR) <- c("w", "FPR");
w_FPR$w <- w_change;

pb <- txtProgressBar(min = 1, max = length(w_change));
for(i in 1:length(w_change)){

	w_ix <- w_change[i];

	# Calculate fpr
	options(warn = -1) ;
	fp <- NULL;
	for(m in 1:10000){
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


## 4.b) Figure 2: TPR vs w, setting TE = 0.7 (JR)
w_TPR <- as.data.frame(matrix(rep(0, length(w_change)*2), nrow = length(w_change), ncol = 2));
names(w_TPR) <- c("w", "TPR");
w_TPR$w <- w_change;

pb <- txtProgressBar(min = 1, max = length(w_change));
for(i in 1:length(w_change)){

	w_ix <- w_change[i];

	# Calculate tpr
	options(warn = -1) ;
	tp <- NULL;
	for(m in 1:10000){
		p_eff <- trial1(TE = 0.7, w = w_ix);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		tp = c(tp, check);

	}; 

	# Incorporate into data.frame
	w_TPR[which(w_TPR$w == w_ix), "TPR"] <- mean(tp, na.rm = T);
	setTxtProgressBar(pb, i);
};


## 4.c) Figure 3: TPR vs TE, setting w = 0.7 (JR)
TPR_TE <- as.data.frame(matrix(rep(0, length(te_change)*2), nrow = length(te_change), ncol = 2));
names(TPR_TE) <- c("TE", "TPR");
TPR_TE$TE <- te_change;


pb <- txtProgressBar(min = 1, max = length(te_change));
for(i in 1:length(te_change)){

	te_ix <- te_change[i];

	# Calculate tpr
	options(warn = -1) ;
	tp <- NULL;
	for(m in 1:10000){
		p_eff <- trial1(TE = te_ix, w = 0.7);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		tp = c(tp, check);

	}; 

	# Incorporate into data.frame
	TPR_TE[which(TPR_TE$TE == te_ix), "TPR"] <- mean(tp, na.rm = T);
	setTxtProgressBar(pb, i);
};


## 4.d) Plotting (JR)
pdf(file = "~/courses/epib635/graphs/epib635_YL_dalGen_mkd_JR_scatterplot_permutations_single_arm_1500.pdf");

# Figure 1
plot <- ggplot(w_FPR, aes(x = w, y = FPR)) + 
	geom_line(size = 1) +
	labs(title = "False positive rate vs Weight",
		subtitle = "Historical control borrowing",
		x = "Weight",
		y = "False positive rate"
		) + 
	scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
	theme_bw() +
	theme(
		panel.border = element_rect(linetype = "solid", colour = "black", size = 2), 
		axis.text = element_text(size = 15, angle = 30, hjust = 1),
		axis.title = element_text(size = 17),
		axis.ticks = element_line(size = 1),
		plot.title = element_text(size = 20),
		plot.subtitle = element_text(size = 15)
		);
print(plot);	

# Figure 2
plot <- ggplot(w_TPR, aes(x = w, y = TPR)) + 
	geom_line(size = 1) +
	labs(title = "Power vs Weight",
		subtitle = "Historical control borrowing",
		x = "Weight",
		y = "Power"
		) + 
	scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue", size = 1) +
	theme_bw() +
	theme(
		panel.border = element_rect(linetype = "solid", colour = "black", size = 2), 
		axis.text = element_text(size = 15, angle = 30, hjust = 1),
		axis.title = element_text(size = 17),
		axis.ticks = element_line(size = 1),
		plot.title = element_text(size = 20),
		plot.subtitle = element_text(size = 15)
		);
print(plot);	

# Figure 3
plot <- ggplot(TPR_TE, aes(x = TE, y = TPR)) + 
	geom_line(size = 1) +
	labs(title = "Power vs Treatment efficacy",
		subtitle = "Historical control borrowing, fixed Weight (0.7)",
		x = "Treatment efficacy",
		y = "Power"
		) + 
	scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.1), limits = c(0.3, 0.9)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue", size = 1) +
	theme_bw() +
	theme(
		panel.border = element_rect(linetype = "solid", colour = "black", size = 2), 
		axis.text = element_text(size = 15, angle = 30, hjust = 1),
		axis.title = element_text(size = 17),
		axis.ticks = element_line(size = 1),
		plot.title = element_text(size = 20),
		plot.subtitle = element_text(size = 15)
		);
print(plot);	


dev.off();



# 5) Borrowing for both control and treatment arms
## 5.a) Updating w_bar function to include treatment data (YL)
w_bar = function(w, y, n, a, b, a0 = 1, b0 = 1) {
	f = beta(a + y, b + n - y)/beta(a, b)
	f0 = beta(a0 + y, b0 + n - y)/beta(a0, b0)
	w_bar = w*f/(w*f + (1 - w)*f0)
	return(w_bar)
}

## 5.b) Updating trial simulation function accordingly (YL)
options(warn = -1)
trial2 = function(TE, w,
                  ni = c(99, 198),
                  # keep the same total no. of cases at interim and final analyses
                  a0 = 1, b0 = 1) {
  n = c(ni[1], diff(ni, 1))
  theta = (1 - TE)/(2 - TE)
  p_eff = NULL
  p_success = NULL
  y_i = 0
  i = 0

  repeat {
    i = i + 1
    y_i = y_i + rbinom(1, n[i], theta) # y_i ~ Binom(n_i, theta)
    ai = 1 + y_i
    bi = 1 + ni[i] - y_i

    ### Add robust MAP prior
    # posterior weight:
    wb = w_bar(w, y_i, n[i], ai, bi)

    # P(TE > 0|data):
    # In power_and_mixture_priors.html: y1 is # of cases in historical trial (ctrl arm ) , and y is # of cases in current trial (ctrl arm )

    peff_new = wb*pbeta(0.5,
                        ai + 226 + 180, # add 180 cases in trx to y1 in MAP code
                        bi + (226 + 180) - 226) +
                    (1-wb)*pbeta(0.5,
                                 ai,
                                 bi)
    p_eff = c(p_eff, peff_new)


    p_success = c(p_success,
                   wb*pbetabinom.ab(q = 88,
                                    size = 198,
                                    shape1 = ai + 226 + 180,
                                    # add 180 cases in trx to y1
                                    shape2 = bi + (226 + 180) - 226,
                                    log = FALSE) +
                     (1-wb)*pbetabinom.ab(q = 88,
                                          size = 198,
                                          shape1 = ai,
                                          shape2 = bi,
                                          log = FALSE)
                   )



    if (i <= 1) {
      # P(TE > 0|data) > 0.995;
      if (peff_new > 0.995 | p_success[length(p_success)] < 0.05) break
    } else break
  }
  out = list(round(p_eff, 2), round(p_success, 8))
  return(out)
}

#trial2(TE = 0.7, w = 0.7)
#trial2(TE = 0.7, w = 0.5) # decrease w because trx arm in historical trial may not be reliable


## 5.b) False Positive Rate (YL)
options(warn = -1)
fp = NULL
p_eff <- p_eff_unlist <- check <- list()
for (m in 1:10000) {
  p_eff = trial2(TE = 0, w = 0.7)
  p_eff_unlist = unlist(p_eff[[1]])
  check = ifelse(length(p_eff_unlist) < 2,
                 p_eff_unlist[length(p_eff_unlist)] > 0.995,
                 p_eff_unlist[length(p_eff_unlist)] > 0.986)
  fp = c(fp, check)
}
#head(fp)
#mean(fp, na.rm = T)
#sum(is.na(fp))



## 5.c) True Positive Rate (YL)
# w = 0.7: power not high enough
options(warn = -1)
tp = NULL
for (m in 1:10000) {
  p_eff = trial2(TE = 0.7, w = 0.7)
  p_eff_unlist = unlist(p_eff[[1]])
  check = ifelse(length(p_eff_unlist) < 2,
                 p_eff_unlist[length(p_eff_unlist)] > 0.995,
                 p_eff_unlist[length(p_eff_unlist)] > 0.986)
  tp = c(tp, check)
}
#mean(tp, na.rm = T)
#sum(is.na(tp))


# w = 0.5: power is 1
options(warn = -1)
tp = NULL
for (m in 1:10000) {
  p_eff = trial2(TE = 0.7, w = 0.5)
  p_eff_unlist = unlist(p_eff[[1]])
  check = ifelse(length(p_eff_unlist) < 2,
                 p_eff_unlist[length(p_eff_unlist)] > 0.995,
                 p_eff_unlist[length(p_eff_unlist)] > 0.986)
  tp = c(tp, check)
}
#mean(tp, na.rm = T)
#sum(is.na(tp))


# 6) Identifying relationship between w, TE, and FPR/TPR using historical data from both arms
w_change_new <- seq(0, 1, 0.01);
te_change_new <- seq(0.3, 0.9, 0.01);

## 6.a) Figure 4: FPR vs w, setting TE = 0.0 (JR)
w_FPR_new <- as.data.frame(matrix(rep(0, length(w_change_new)*2), nrow = length(w_change_new), ncol = 2));
names(w_FPR_new) <- c("w", "FPR");
w_FPR_new$w <- w_change_new;

pb <- txtProgressBar(min = 1, max = length(w_change_new));
for(i in 1:length(w_change_new)){

	w_ix <- w_change_new[i];

	# Calculate fpr
	options(warn = -1) ;
	fp <- NULL;
	for(m in 1:10000){
		p_eff <- trial2(TE = 0, w = w_ix);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		fp = c(fp, check);

	}; 

	# Incorporate into data.frame
	w_FPR_new[which(w_FPR_new$w == w_ix), "FPR"] <- mean(fp, na.rm = T);
	setTxtProgressBar(pb, i);
};



## 6.b) Figure 5: TPR vs w, setting TE = 0.7 (JR)
w_TPR_new <- as.data.frame(matrix(rep(0, length(w_change_new)*2), nrow = length(w_change_new), ncol = 2));
names(w_TPR_new) <- c("w", "TPR");
w_TPR_new$w <- w_change_new;

pb <- txtProgressBar(min = 1, max = length(w_change_new));
for(i in 1:length(w_change_new)){

	w_ix <- w_change_new[i];

	# Calculate tpr
	options(warn = -1) ;
	tp <- NULL;
	p_eff <- p_eff_unlist <- check <- list();
	for(m in 1:10000){
		p_eff <- trial1(TE = 0.7, w = w_ix);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		tp = c(tp, check);

	}; 

	# Incorporate into data.frame
	w_TPR_new[which(w_TPR$w == w_ix), "TPR"] <- mean(tp, na.rm = T);
	setTxtProgressBar(pb, i);
};


## 6.c) Figure 6: TPR vs TE, setting w = 0.7 (JR)
TPR_TE_new <- as.data.frame(matrix(rep(0, length(te_change_new)*2), nrow = length(te_change_new), ncol = 2));
names(TPR_TE_new) <- c("TE", "TPR");
TPR_TE_new$TE <- te_change_new;


pb <- txtProgressBar(min = 1, max = length(te_change_new));
for(i in 1:length(te_change_new)){

	te_ix <- te_change_new[i];

	# Calculate tpr
	options(warn = -1) ;
	tp <- NULL;
	for(m in 1:10000){
		p_eff <- trial1(TE = te_ix, w = 0.7);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		tp = c(tp, check);

	}; 

	# Incorporate into data.frame
	TPR_TE_new[which(TPR_TE_new$TE == te_ix), "TPR"] <- mean(tp, na.rm = T);
	setTxtProgressBar(pb, i);
};

## 6.d) Figure 7: TPR vs TE, setting w = 0.5 (JR)
TPR_TE_new.2 <- as.data.frame(matrix(rep(0, length(te_change_new)*2), nrow = length(te_change_new), ncol = 2));
names(TPR_TE_new.2) <- c("TE", "TPR");
TPR_TE_new.2$TE <- te_change_new;


pb <- txtProgressBar(min = 1, max = length(te_change_new));
for(i in 1:length(te_change_new)){

	te_ix <- te_change_new[i];

	# Calculate tpr
	options(warn = -1) ;
	tp <- NULL;
	p_eff <- p_eff_unlist <- check <- list();
	for(m in 1:10000){
		p_eff <- trial1(TE = te_ix, w = 0.7);
		p_eff_unlist  <- unlist(p_eff[[1]]);
		check <- ifelse(length(p_eff_unlist) < 2,
			 p_eff_unlist[length(p_eff_unlist)] > 0.995,
			 p_eff_unlist[length(p_eff_unlist)] > 0.986);
		tp = c(tp, check);

	}; 

	# Incorporate into data.frame
	TPR_TE_new.2[which(TPR_TE_new.2$TE == te_ix), "TPR"] <- mean(tp, na.rm = T);
	setTxtProgressBar(pb, i);
};


## 6.d) Plotting (JR)
pdf(file = "~/courses/epib635/graphs/epib635_YL_dalGen_mkd_JR_scatterplot_permutations_both_arms_1500.pdf");

# Figure 4
plot <- ggplot(w_FPR_new, aes(x = w, y = FPR)) + 
	geom_line(size = 1) +
	labs(title = "False positive rate vs Weight",
		subtitle = "Historical control + treatment borrowing",
		x = "Weight",
		y = "False positive rate"
		) + 
	scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
	theme_bw() +
	theme(
		panel.border = element_rect(linetype = "solid", colour = "black", size = 2), 
		axis.text = element_text(size = 15, angle = 30, hjust = 1),
		axis.title = element_text(size = 17),
		axis.ticks = element_line(size = 1),
		plot.title = element_text(size = 20),
		plot.subtitle = element_text(size = 15)
		);
print(plot);	

# Figure 5
plot <- ggplot(w_TPR_new, aes(x = w, y = TPR)) + 
	geom_line(size = 1) +
	labs(title = "Power vs Weight",
		subtitle = "Historical control + treatment borrowing",
		x = "Weight",
		y = "Power"
		) + 
	scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue", size = 1) +
	theme_bw() +
	theme(
		panel.border = element_rect(linetype = "solid", colour = "black", size = 2), 
		axis.text = element_text(size = 15, angle = 30, hjust = 1),
		axis.title = element_text(size = 17),
		axis.ticks = element_line(size = 1),
		plot.title = element_text(size = 20),
		plot.subtitle = element_text(size = 15)
		);
print(plot);	

# Figure 6
plot <- ggplot(TPR_TE_new, aes(x = TE, y = TPR)) + 
	geom_line(size = 1) +
	labs(title = "Power vs Treatment efficacy",
		subtitle = "Historical control + treatment borrowing, fixed Weight (0.7)",
		x = "Treatment efficacy",
		y = "Power"
		) + 
	scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.1), limits = c(0.3, 0.9)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue", size = 1) +
	theme_bw() +
	theme(
		panel.border = element_rect(linetype = "solid", colour = "black", size = 2), 
		axis.text = element_text(size = 15, angle = 30, hjust = 1),
		axis.title = element_text(size = 17),
		axis.ticks = element_line(size = 1),
		plot.title = element_text(size = 20),
		plot.subtitle = element_text(size = 15)
		);
print(plot);	


# Figure 7
plot <- ggplot(TPR_TE_new.2, aes(x = TE, y = TPR)) + 
	geom_line(size = 1) +
	labs(title = "Power vs Treatment efficacy",
		subtitle = "Historical control + treatment borrowing, fixed Weight (0.5)",
		x = "Treatment efficacy",
		y = "Power"
		) + 
	scale_x_continuous(breaks = seq(0.3, 0.9, by = 0.1), limits = c(0.3, 0.9)) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
	geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue", size = 1) +
	theme_bw() +
	theme(
		panel.border = element_rect(linetype = "solid", colour = "black", size = 2), 
		axis.text = element_text(size = 15, angle = 30, hjust = 1),
		axis.title = element_text(size = 17),
		axis.ticks = element_line(size = 1),
		plot.title = element_text(size = 20),
		plot.subtitle = element_text(size = 15)
		);
print(plot);	
dev.off();

stop.time <- Sys.time();
end.time <- stop.time-start.time
# For 10,000 simulations:
# Time difference of 2.735275 hours

# 7) Save image
savey <- c("~/courses/epib635/rdata_images/epib635_YL_dalGen_mkd_JR_modifications_1500.RData");
load(savey);

save.image(savey);

