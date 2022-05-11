fst = function(project,run = 1, K, ploidy = 2){
library(LEA)
l = dim(G(project, K = K, run = run))[1]
q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
if (ploidy == 2) {
G1.t = G(project, K = K, run = run)[seq(2,l,by = 3),]
G2.t = G(project, K = K, run = run)[seq(3,l,by = 3),]
freq = G1.t/2 + G2.t}
else {
freq = G(project, K = K, run = run)[seq(2,l,by = 2),]}
H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
H.t = P.t*(1-P.t)
return(1-H.s/H.t)
}
