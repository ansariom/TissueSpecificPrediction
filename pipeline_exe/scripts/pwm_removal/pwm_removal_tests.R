
df <- data.frame(a = c(), coef = c(), type = c())
dir <- "~/Downloads/ibdc/coefs/"

top0 <- read.table("~/Downloads/ibdc/coefs/top0_coefs.txt", col.names = c("feature", "coef"))
top0$rank0 <- as.numeric(row.names(top0))
top10 <- read.table("~/Downloads/ibdc/coefs/top10_coefs.txt", col.names = c("feature", "coef"))
top10$rank10 <- as.numeric(row.names(top10))
top <- merge(top0, top10, by = "feature")
top$rank_diff <- top$rank10 - top$rank0
top <- top[order(as.numeric(top$rank0)),]
top <- top[1:200,]
hist(top$rank_diff)


files <- list.files(dir)
i = 1
x <- c(0, 10, 20, 50, 100, 150)

high_diff <- data.frame(feature=c(), diff=c()) 
for (i in 1:(length(x)-1)) {
  for (j in seq(i+1,length(x))) {
    #print(paste(i, j , sep = "-"))
    f1 = paste(dir,"/top", x[i], "_coefs.txt", sep = "")
    f2 = paste(dir, "/top", x[j], "_coefs.txt", sep = "")
    print(f1)
    print(f2)
    a <- read.table(f1, header = F, col.names = c("feature", "coef"))
    b <- read.table(f2, header = F, col.names = c("feature", "coef"))
    a$rank1 <- as.numeric(rownames(a))
    b$rank2 <- as.numeric(rownames(b))
    m <- merge(a, b, by = "feature")
    m <- m[order(m$rank1),]
    m <- m[1:300,]
    m$diff <- abs(m$rank1 - m$rank2)
    high_diff <- rbind(high_diff, m[m$diff > 5000, c("feature", "diff")])
    f <- paste(dir, "/top", x[i], "_top", x[j], ".jpg", sep = "")
    jpeg(filename = f)
    plot(m$rank1, m$rank2, xlab = "rank before PWM removal", ylab = "rank after PWM removal", 
         main = paste("top", x[i], " and top", x[j], sep = ""))
    dev.off()
  }
}


ranked_df <- data.frame(feature=c(), rank = c(), type = c())

i = 1
for (f in files) {
  print(f)
  n <- x[i]
  print(n)
  i <- i +1
  s <- paste(dir, "/", f, sep = "")
  a <- read.table(s, header = F, col.names = c("feature", "coef"))
  a$type <- paste("top", n, sep = "")
  a$rank <- row.names(a)
  a$coef <- abs(a$coef)
  df <- rbind(df, a[1:500, c("rank", "coef", "type")])
  ranked_df <- rbind(ranked_df, a[1:10, c("feature", "rank", "type")])
  
}

library(ggplot2)
ggplot(df, aes(as.numeric(rank), coef)) + geom_point() + facet_wrap(~type)
ggplot(ranked_df, aes(feature,as.numeric(rank))) + geom_point() + facet_wrap(~type)
