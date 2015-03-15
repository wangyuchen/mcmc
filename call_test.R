# subroutine mcmc(l,g0,g1,g,p,result)
# integer, intent(in)::l, g0, g1, g
# real(kind = 8), dimension(0:2, l), intent(in) :: p
# integer, dimension(l, 0:2), intent(out):: result

setwd("~/R/mcmc/")

p <- read.table("p_matrix.txt")
p <- t(p)
is.double(p)

mcmc <- function(p, l = 2000L, g0 = 100L, g1 = 200L, g = l * 1000L) {
  res <- matrix(0L, l, 3)
  dyn.load("mcmc_sub.so")
  sum_counts <- .Fortran("mcmc", l, g0, g1, g, p, res)[[6]]
  return(sum_counts)
}

sum_counts <- mcmc(p)
Iprop <- apply(sum_counts, MARGIN = 1, which.max) - 1
summary(as.factor(Iprop))


result_Iprop3 <- read.table("~/R/mcmc/result_Iprop3.txt", quote="\"")
result_Iprop3 <- as.vector(t(result_Iprop3))
summary(as.factor(result_Iprop3))

trueI <- read.table("./trueI.txt")
trueI <- trueI$V1
summary(as.factor(trueI))

result_Iprop <- scan()
summary(as.factor(result_Iprop))


sum(trueI == result_Iprop)
sum(trueI == Iprop)
sum(result_Iprop == Iprop)
