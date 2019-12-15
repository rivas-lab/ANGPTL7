b<- matrix(nrow =  810 , ncol = 2)
e <- matrix(nrow = 155656, ncol = 2)
d<- matrix(nrow = 13810 , ncol = 2)
a<- matrix(nrow = 86, ncol = 2)
e[,] <- 0
d[,1] <- 0
d[,2] <- 1
b[,1] <- 1
b[,2] <- 0
a[,] <- 1
df.1 <- rbind(a,b,d,e)
a <- c(rep(1,2),rep(0,86-2))
b <- c(rep(1,57),rep(0,810-57))
d <- c(rep(1,183),rep(0,13810 - 183))
e <- c(rep(1,3133),rep(0,155656 - 3133))
df.2 <- c(a,b,d,e)
colnames(df.2) <- c("myoc","angptl7")
colnames(df.1) <- c("myoc","angptl7")
summary(glm(df.2 ~ df.1[,"myoc"]*df.1[,"angptl7"], family = binomial))
savehistory(file = "angptl7myoc.interaction.txt")
