#############################################################
"
This is a script to run the power analysis on the association of germline 
mutations and somatic mutations in certain frequencies.

Question: I have a germline mutation in frequency 1% and a somatic mutation 
in frequency 10%. Whats the power I have across different sizes of cohort. 

Note that we already know that the germline and the somatic mutation are present 
in the cohort and we don't model the probabilities to find these mutations in a 
cohort of certain size.

"
#############################################################

## We do it via simulation
res = NULL
ps = list(c(0.1, 0.01), c(0.1, 0.1), c(0.2, 0.15), c(0.3, 0.20))
effect_sizes = c(1.2, 1.5, 1.8, 2, 2.5)
for(p in ps){
    print(p)
    for(e in effect_sizes){
        print(e)
        for(y in 100:1000)
        {   
            bootTest <- rep(FALSE, 1000) #store results of thousand simulations
            n <- y  # guess sample size - we could iterate/search across values, I just changed until I got roughly 80% power
            p1 <- p[2]   # frequency of mutations
            p2 <- p[1]  # frequency of 'associated' mutations
            effectSize <- e #hypothesised
            for (i in seq(along=bootTest)) { #run many simulations
                r1 <- runif(n)<p1 # randomly generate TRUEs with correct frequency
                if(sum(r1)==0){ ## If there is no TRUE in the vector test gives an error - add a random TRUE
                    r1[sample(1:n, 1)] = TRUE
                }
                r2 <- runif(n)<ifelse(r1, p2*effectSize, p2) # create correlated TRUEs
                if(sum(r2)==0){ ## If there is no TRUE in the vector test gives an error - add a random TRUE
                    r2[sample(1:n, 1)] = TRUE
                }
                bootTest[i] <- chisq.test(table(r1, r2))$p.value < 0.05 # store whether true-positive is found
            }
            m =mean(bootTest) # will be the proportion of true positives in total positives = power
            res = rbind(res, data.frame(p1=p1, p2=p2, effect_size=e, power=m, n=y))
        }
    }
}

res2 = res2 %>% mutate(grp=paste0(p1, "_", p2, "_", effect_size))
t = res1

library(RColorBrewer)
n <- length(unique(t$grp))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(t, aes(x=n, y=power, group=grp, color=grp)) + 
    geom_smooth() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    scale_x_log10(breaks=c(seq(100, 1000, 100))) +
    theme_boss() +
    theme(
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")
    ) +
    scale_color_manual(values = col_vector)



library(powerGWASinteraction)

res = NULL
for (y in 10:1000){
    print(y)
    mod1 <- list(prev=1, pGene1=0.1, pGene2=0.1, beta.LOR=c(0,0,1), nSNPs=3)
    r1 = powerGG(n=y,mod=mod1,caco=0.01,alpha=.05, alpha1 = 1)[["case.control"]]
    res = rbind(res, data.frame(p1=0.5, p2=0.5, power=r1, n=y))
}
res %>% head








