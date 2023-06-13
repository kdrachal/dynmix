

### object - mixest obtained from mixest1()
### x.post - matrix T x m of independent variables for post-intervention period
### m      - number of independent variables
### y.post - vector of post-intervention period dependent variable
### alpha  - desired tail area probability for posterior intervals
### n.sim  - number of post-intervention period simulations


cauimp <- function(object,x.post,y.post,alpha=0.05,n.sim=100)
  {
    T <- length(y.post)
    
    y.sim <- matrix(0,ncol=n.sim,nrow=T)
    for (i in 1:n.sim)
      {
        y.sim[,i] <- .mixest1.sim(object=object,x=x.post)
      }
    
    point.pred <- rowMeans(y.sim)
    
    p.down <- alpha / 2
    p.up <- 1 - (alpha / 2)
    
    y.post.rep <- matrix(y.post,ncol=n.sim,nrow=T,byrow=FALSE)    
    
    rel.eff <- (colSums(y.post.rep) / colSums(y.sim)) - 1
    
    summary <- matrix(NA,ncol=2,nrow=13)
    colnames(summary) <- c("average","cumulative")
    rownames(summary) <- c("observed","predicted","predicted U","predicted D","predicted sd",
                           "absolute effect","absolute effect U","absolute effect D","absolute effect sd",
                           "relative effect","relative effect U","relative effect D","relative effect sd")
    summary[1,1] <- mean(y.post)
    summary[1,2] <- sum(y.post)
    summary[2,1] <- mean(point.pred)
    summary[2,2] <- sum(point.pred)
    summary[3,1] <- quantile(colMeans(y.sim),p.up)
    summary[3,2] <- quantile(colSums(y.sim),p.up)
    summary[4,1] <- quantile(colMeans(y.sim),p.down)
    summary[4,2] <- quantile(colSums(y.sim),p.down)
    summary[5,1] <- sd(colMeans(y.sim))
    summary[5,2] <- sd(colSums(y.sim))
    summary[6,1] <- summary[1,1] - summary[2,1]
    summary[6,2] <- summary[1,2] - summary[2,2]
    summary[7,1] <- quantile(colMeans(y.post.rep - y.sim),p.up)
    summary[7,2] <- quantile(colSums(y.post.rep - y.sim),p.up)
    summary[8,1] <- quantile(colMeans(y.post.rep - y.sim),p.down)
    summary[8,2] <- quantile(colSums(y.post.rep - y.sim),p.down)
    summary[9,1] <- sd(colMeans(y.post.rep - y.sim))
    summary[9,2] <- sd(colSums(y.post.rep - y.sim))
    summary[10,1] <- mean(rel.eff) 
    summary[10,2] <- summary[10,1]
    summary[11,1] <- quantile(rel.eff,p.up)
    summary[11,2] <- summary[11,1]
    summary[12,1] <- quantile(rel.eff,p.down)
    summary[12,2] <- summary[12,1]
    summary[13,1] <- sd(rel.eff)
    summary[13,2] <- summary[13,1]

    p <- min(sum(c(colSums(y.sim),summary[1,2]) >= summary[1,2]),
             sum(c(colSums(y.sim),summary[1,2]) <= summary[1,2]))
    p <- p / (n.sim + 1)
    
    sig <- (!((summary[12,1] < 0) & (summary[11,1] > 0))) 

    out <- list(summary,sig,p,point.pred,alpha,n.sim,y.sim)
    names(out) <- c("statistics","significance","p","y.hat","alpha","n.sim","y.sim")
    return(out)
  }
   
   