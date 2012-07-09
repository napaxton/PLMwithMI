m1 <- pdata.frame(m1, "country", time="year")
m2 <- pdata.frame(m2, "country", time="year")
m3 <- pdata.frame(m3, "country", time="year")
m4 <- pdata.frame(m4, "country", time="year")
m5 <- pdata.frame(m5, "country", time="year")

MIlist <- list(m1,m2,m3,m4,m5)

save(MIlist,file="~/Desktop/MIlist.Rdata")

plm.ran<-list()						# List of plm objects
plm.with<-list()					# "
results.random<-list()				# List of results objects
results.within<-list()				# "
summary.plm.out.ran <-list()		# List of summary objects
summary.plm.out.with <-list()		# "
effects.random <-list()				#
summary.ran.tss <-list()			# Total sum of squares objects
summary.ran.ssr <-list()			# Sum of squared residuals
summary.ran.F <-list()				# F-test
summary.random.resduals <- list()	# Residuals
summary.with.tss <-list()
summary.with.ssr <-list()
summary.with.F <-list()
summary.within.residuals <- list()
results.FE <- list()
results.hausman <-list()			# Hausman test results
summary.random.residuals <- list()	# Summary of the residuals

form <- y ~ x1 + x2 +x3

for (i in 1:length(MIlist)){#this loop does the model for each of the datasets in the list MIlist.
	plm.ran[[i]] <- plm(form, data=MIlist[[i]],model="random")
	plm.with[[i]] <- plm(form, data=MIlist[[i]],model="within")
	results.random[[i]] <- plm.ran[[i]] #Model and coefs of the random effects model
	results.within[[i]] <- plm.with[[i]]

	summary.plm.out.ran[[i]] <- summary(plm.ran[[i]]) #Summary of the effects, residuals, coefs and overall stats of the random effects model
	summary.plm.out.with[[i]] <- summary(plm.with[[i]]) #Summary of the coefs and overall test stats of the within effects model

	summary.random.residuals[[i]] <- summary(plm.ran[[i]])$residuals
	effects.random[[i]] <- summary(plm.ran[[i]])$sigma2 #Model and coefs of the random effects model
	summary.ran.tss[[i]] <- summary(plm.ran[[i]])$tss #Total sum squares of the random model
	summary.ran.ssr[[i]] <- summary(plm.ran[[i]])$ssr #Total sum squares of the random model
	summary.ran.F[[i]] <- summary(plm.ran[[i]])$fstatistic$statistic #Total sum squares of the random model

	summary.within.residuals[[i]] <- summary(plm.with[[i]])$residuals
	summary.with.tss[[i]] <- summary(plm.with[[i]])$tss #Total sum squares of the random model
	summary.with.ssr[[i]] <- summary(plm.with[[i]])$ssr #Total sum squares of the random model
	summary.with.F[[i]] <- summary(plm.with[[i]])$fstatistic$statistic #Total sum squares of the random model	results.within[[i]] <- plm.out[[i]]$within

	# results.FE[[i]] <- FE(plm.out[[i]])
	results.hausman[[i]] <- phtest(plm.with[[i]],plm.ran[[i]]) #Hausman test results
}

#We need the degrees of freedom for each model for the F stats hypo testing later, and since they are the same for each run, they can be acquired and stored outside of the above loop.
ran.df1 <- summary(plm.ran[[i]])$fstatistic$parameter[1]
ran.df2 <- summary(plm.ran[[i]])$fstatistic$parameter[2]
with.df1 <- summary(plm.with[[i]])$fstatistic$parameter[1]
with.df2 <- summary(plm.with[[i]])$fstatistic$parameter[2]

# *** STORE THE COEFFICIENTS OF THE REGRESSIONS ON THE IMPUTATIONS IN THEIR OWN MATRICES ***
coef.store.ran <- matrix(NA, ncol=length(results.random[[i]]$coef), nrow=length(MIlist))#container for coefficients
coef.store.with <- matrix(NA, ncol=length(results.within[[i]]$coef), nrow=length(MIlist))#container for coefficients

for (i in 1:length(results.random)){#gathering the coefficients from each of the results
	coef.store.ran[i,] <- results.random[[i]]$coef #putting all of the coefs in one matrix.
}

for (i in 1:length(results.within)){#gathering the coefficients from each of the results
	coef.store.with[i,] <- results.within[[i]]$coef #putting all of the coefs in one matrix.
}

# *** COEFFICIENT POINT ESTIMATES (MEANS OF THE MATRIX OF IMPUTATION STORAGE MATRICES)
coef.final.ran <-apply(coef.store.ran, 2, mean)#getting the point estimates by taking the means of the coefficients
names(coef.final.ran) <- names(results.random[[1]]$coefficients)
coef.final.ran

coef.final.with <-apply(coef.store.with, 2, mean)#getting the point estimates by taking the means of the coefficients
names(coef.final.with) <- names(results.within[[1]]$coefficients)
coef.final.with

# *** STORE THE Standard errors OF THE REGRESSIONS ON THE IMPUTATIONS IN THEIR OWN MATRICES ***
se.store.ran<-matrix(NA, ncol=length(diag(results.random[[i]]$vcov)), nrow=length(MIlist)) #container for the std errs.
for (i in 1:length(results.random)){
	se.store.ran[i,]<-diag(vcov(results.random[[i]])) #storing the variances in the container
}

se.store.with<-matrix(NA, ncol=length(diag(results.within[[i]]$vcov)), nrow=length(MIlist)) #container for the std errs.
for (i in 1:length(results.within)){
	se.store.with[i,]<-diag(vcov(results.within[[i]])) #storing the variances in the container
}

se.1.ran<-apply(se.store.ran, 2, mean)#taking the mean of the variances.
se.1.with<-apply(se.store.with, 2, mean)#taking the mean of the variances.
se.2.ran<-apply(se.store.ran, 2, var)*(1+(1/length(results.random)))#taking the variance of the variance point estimates from the mean
se.2.with<-apply(se.store.with, 2, var)*(1+(1/length(results.within)))#taking the variance of the variance point estimates from the mean
se.final.ran<-sqrt(se.1.ran+se.2.ran)#calculating the final standard erros.
se.final.with<-sqrt(se.1.with+se.2.with)#calculating the final standard erros.

t.stat.ran <- coef.final.ran/se.final.ran
t.stat.with <- coef.final.with/se.final.with

pr.gt.t.ran <- dt(t.stat.ran, results.random[[i]]$df.residual)
pr.gt.t.with <- dt(t.stat.with, results.within[[i]]$df.residual)

# OVERALL STATISTICS 
effects.store.random <- matrix(NA, ncol=(length(effects.random[[i]])-1), nrow=length(MIlist)) #container for the model effects statistics
for (i in 1:length(effects.random)){
	effects.store.random[i,1] <- effects.random[[i]]$idios
	effects.store.random[i,2] <- effects.random[[i]]$id
}
effects.1.ran <-apply(effects.store.random, 2, mean) #have to do in two steps like the SE above, because these are variances
effects.2.ran <- apply(effects.store.random,2,var)*(1+(1/length(results.random)))
effects.ran.var <-effects.1.ran+effects.2.ran
effects.ran.sd <-sqrt(effects.ran.var)
effects.ran.share <- effects.ran.var/sum(effects.ran.var)

overall.stats.random<- matrix(NA, nrow=3,ncol=length(summary.plm.out.ran))

# This loop extracts the random model's anova statistics
for (i in 1:length(summary.plm.out.ran)){
	overall.stats.random[1,i] <- summary.ran.tss[[i]]
	overall.stats.random[2,i] <- summary.ran.ssr[[i]]
	overall.stats.random[3,i] <- summary.ran.F[[i]]
}
overall.ran.tss <-mean(overall.stats.random[1,])
overall.ran.ssr <-mean(overall.stats.random[2,])
overall.ran.F <-mean(overall.stats.random[3,])
overall.ran.sig <- df(overall.ran.F,ran.df1,ran.df2)

# This loop extracts the within model's anova statistics
overall.stats.within<- matrix(NA, nrow=3,ncol=length(summary.plm.out.with))
for (i in 1:length(summary.plm.out.with)){
	overall.stats.within[1,i] <- summary.with.tss[[i]]
	overall.stats.within[2,i] <- summary.with.ssr[[i]]
	overall.stats.within[3,i] <- summary.with.F[[i]]
}
overall.with.tss <-mean(overall.stats.within[1,])
overall.with.ssr <-mean(overall.stats.within[2,])
overall.with.F <-mean(overall.stats.within[3,])
overall.with.sig <- df(overall.with.F,with.df1,with.df2)

# *** CREATE A MATRIX TO STORE THE CHI^2 STATITICS FOR EACH IMPUTATION'S HAUSMAN TEST SCORE
hausman.store<- matrix(NA, ncol=1, nrow=length(results.hausman))
for (i in 1:length(results.hausman)){
	hausman.store[i,] <- results.hausman[[i]]$statistic #cbind(c("Hausman for imputation ",i),
}

hausman.final <- apply(hausman.store,2,mean) #taking the means of the chi^2 statistics
hausman.report<- matrix(NA, nrow=length(results.hausman) + 1, ncol=1)
for (i in 1:length(results.hausman)){
	hausman.report[i,] <- results.hausman[[i]]$statistic
}

hausman.report[length(results.hausman)+1,]<-hausman.final #CREATE A REPORT MATRIX FOR THE SET OF HAUSMAN TESTS (SINCE I WANT TO SHOW THE INDIVIDUAL RESULTS, AS WELL AS THE AVERAGED SCORE)
colnames(hausman.report) <- cat(paste("Chi-sqared with ",results.hausman[[i]]$parameter," df"))
data.frame(hausman.report)
hausmanrows<-matrix(NA, nrow=length(hausman.report)-1, ncol=1)
for (i in 1:length(hausman.report)-1){
	hausmanrows[i,]<- i
}
hausmanrows <- append(hausmanrows, "Average of all imputations") #, after=5)
rownames(hausman.report) <- hausmanrows

# *** COMBINE AND OUTPUT THE FINAL RESULTS ***

sink(file="~/Desktop/plmresults.txt", split=TRUE)

#RANDOM EFFECTS MODEL
cat("========== RANDOM EFFECTS MODEL SUMMARY ==========\n")
cat("\nCall:\n")
  print(results.random[[i]]$call)
pdim(plm.ran[[i]])

effects.final <- cbind(var=effects.ran.var, std.dev=effects.ran.sd,share=effects.ran.share)
row.names(effects.final) <- c("Idiosyncratic ", "Individual ")
cat("\nEffects:\n")
printCoefmat(effects.final)
final.results.ran <- cbind(coef.final.ran, se.final.ran, t.stat.ran, pr.gt.t.ran) #putting all the regression results together and printing them out.
row.names(final.results.ran) <- names(results.random[[1]]$coefficients)
colnames(final.results.ran) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
printCoefmat(final.results.ran)
cat(paste("Total Sum of Squares: ",signif(overall.ran.tss,getOption("digits")),"\n",sep=""))
cat(paste("Residual sum of squares: ",signif(overall.ran.ssr,getOption("digits")),"\n",sep=""))
ran.rsq <-  1-(overall.ran.ssr/overall.ran.tss)
cat(paste("R-Squared: ",signif(ran.rsq,getOption("digits")),"\n",sep=""))

cat(paste("F-statistic: ",signif(overall.ran.F),
			" on ",ran.df1," and ",ran.df2,
			" DF, p-value: ",signif(overall.ran.sig),"\n",sep=""))

#WITHIN EFFECTS MODEL
cat("\n\n========== WITHIN EFFECTS MODEL SUMMARY ==========\n")
cat("\nCall:\n")
  print(results.within[[i]]$call)
pdim(plm.with[[i]])

final.results.with <- cbind(coef.final.with, se.final.with, t.stat.with, pr.gt.t.with) #putting all the results together and printing them out.
row.names(final.results.with) <- names(results.within[[1]]$coefficients)
colnames(final.results.with) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
printCoefmat(final.results.with)
cat(paste("Total Sum of Squares: ",signif(overall.with.tss,getOption("digits")),"\n",sep=""))
cat(paste("Residual sum of squares: ",signif(overall.with.ssr,getOption("digits")),"\n",sep=""))
with.rsq <-  1-(overall.with.ssr/overall.with.tss)
cat(paste("R-Squared: ",signif(with.rsq,getOption("digits")),"\n",sep=""))

cat(paste("F-statistic: ",signif(overall.with.F),
			" on ",with.df1," and ",with.df2,
			" DF, p-value: ",signif(overall.with.sig),"\n",sep=""))

#Which to choose, the random or within?
cat("\n\n========== HAUSMAN TEST (WITHIN & RANDOM EFFECTS MODELS) REPORT ==========\n")
hausman.report
cat(paste("p-value: ", signif(dchisq(hausman.final, df=results.hausman[[i]]$parameter))))

# *** COMBINE AND OUTPUT THE FINAL RESULTS *** SECOND WAY

#RANDOM EFFECTS MODEL
results.random$call
pdim(results.random)
effects.final <- cbind(effects.ran.var, effects.ran.sd,effects.ran.share)
row.names(effects.final) <- c("Idiosyncratic ", "Individual ")
effects.final
final.results.ran <- cbind(coef.final.ran, se.final.ran, t.stat.ran, pr.gt.t.ran) #putting all the regression results together and printing them out.
row.names(final.results.ran) <- names(results.random[[1]]$coefficients)
final.results.ran
ran.tss <- c(print("Total sum of squares"), overall.ran.tss)
ran.tss
ran.ssr <- c(print("Residual sum of squares"), overall.ran.ssr)
ran.ssr
ran.rsq <- c(print("Rsq"), (1-overall.ran.ssr/overall.ran.tss))
ran.rsq
ran.fstat <- c(print("F"), overall.ran.F)
ran.fprob <- c(print("P(F>0)"), overall.ran.sig)

#WITHIN EFFECTS MODEL
results.within$call
pdim(results.within)

final.results.with <- cbind(coef.final.with, se.final.with, t.stat.with, pr.gt.t.with) #putting all the results together and printing them out.
row.names(final.results.with) <- names(results.within[[1]]$coefficients)
final.results.with
with.tss <- c(print("Total sum of squares"), overall.with.tss)
with.tss
with.ssr <- c(print("Residual sum of squares"), overall.with.ssr)
with.ssr
with.rsq <- c(print("Rsq"), (1-overall.with.ssr/overall.with.tss))
with.rsq
with.fstat <- c(print("F"), overall.with.F)
with.fprob <- c(print("P(F>0)"), overall.with.sig)
#Which to choose, the random or within?
hausman.report

sink(file="~/Desktop/plmresults.txt", split=TRUE)
Sys.time(); Sys.timezone()
print("These are the results of panel data analysis on five multiply imputed datasets.")
print("===============================================================================")
print("-------------------------------------------------------------------------------")
print("Here's what the summary of the random models looks like:")
print("The panel data random effects regression equation "); print(form)
print("-------------------------------------------------------------------------------")
print("The panel dimensions: "); print(pdim(data1))
print("-------------------------------------------------------------------------------")
print("Overall model effects: "); print(effects.final)
print("-------------------------------------------------------------------------------")
print("Table of coefficients: "); print(final.results.ran)
print("-------------------------------------------------------------------------------")
print("Overall model statistics:")
ran.tss
ran.ssr
ran.rsq
ran.fstat
ran.fprob

print("")
print("")
print("")
print("===============================================================================")
print("Here's what the summary of the within models looks like:")
print("The panel data within effects regression equation "); print(form)
print("-------------------------------------------------------------------------------")
print("The panel dimensions: "); print(pdim(data1))
print("-------------------------------------------------------------------------------")
print("Table of coefficients: "); print(final.results.with)
print("-------------------------------------------------------------------------------")
print("Overall model statistics:")
with.tss
with.ssr
with.rsq
with.fstat
with.fprob
print("")
print("")
print("")
print("===============================================================================")

print("HAUSMAN TEST RESULTS FOR INDIVIDUAL IMPUTATIONS AND THE AVERAGE")

hausman.report
print("")
print("===============================================================================")
sink()
