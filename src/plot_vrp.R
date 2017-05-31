setwd("F:/Universidad/MASTER/algoevo/plot_crossover")

nodes <- read.table("a-n80-k10.vrp")

# Plot all nodes

plot(nodes[1,2], nodes[1,3], xlim = c(min(nodes[,2]), max(nodes[,2])), ylim = c(min(nodes[,3]), max(nodes[,3])), ylab = "Y-Cartesian points", xlab = "X-Cartesian points", main = "Vehicle Routing Problem")
for(i in 1:length(nodes[,1])){
  points(nodes[i,2], nodes[i,3])
}

# plot depot in color
points(nodes[1,2], nodes[1,3], col = 'purple', lwd = 4)



# Plot solutions
rand1 <- c(45, 22, 50, 76, 72, 66, 67, 70, 0, 36, 54, 9, 53, 55, 33, 15, 64, 42, 0, 56, 69, 77, 51, 41, 65, 39, 47, 46, 25, 60, 0, 3, 35, 19, 26, 31, 74, 17, 75, 0, 20, 57, 29, 13, 27, 59, 0, 61, 16, 5, 43, 78, 68, 44, 12, 30, 0, 23, 6, 8, 24, 62, 37, 0, 2, 34, 11, 63, 1, 52, 0, 79, 28, 7, 71, 48, 0, 18, 10, 14, 21, 40, 49, 4, 0, 38, 58, 32, 73)
rand1 <- c(0, rand1)
rand1 <- c(rand1, 0)
#rand2 <- c(28, 40, 33, 2, 21, 15, 0, 7, 37, 8, 39, 14, 32, 45, 30, 43, 17, 6, 27, 5, 36, 18, 26, 16, 42, 29, 35, 19, 11, 10, 22, 13, 24, 12, 46, 20, 31, 38, 47, 4, 41, 23, 9, 44, 34, 3, 25, 1)
rand1 <- rand1+1
#rand2 <- rand2+1

mycolors <- c('red', 'blue', 'yellow', 'black', 'purple', 'green', 'cyan', 'brown', 'orange', 'chartreuse','red', 'blue', 'yellow', 'black', 'purple', 'green', 'cyan', 'brown', 'orange', 'chartreuse')

currcolor <- 1

for(i in 2:length(rand1)){
  
  lines(c(nodes[rand1[i-1],2], nodes[rand1[i],2]) , c(nodes[rand1[i-1],3], nodes[rand1[i],3]), col = mycolors[currcolor])
  if(rand1[i] == 1){
    currcolor <- (currcolor + 1) %% length(mycolors)
  }
  
}
lines(c(nodes[rand1[1],2], nodes[rand1[length(rand1)],2]) , c(nodes[rand1[1],3], nodes[rand1[length(rand1)],3]), col = mycolors[currcolor])

# for(i in 2:length(rand2)){
#   points(c(nodes[rand2[i-1],2], nodes[rand2[i],2]) , c(nodes[rand2[i-1],3], nodes[rand2[i],3]), type = "b", col = "blue", lty="dashed")
# }
# points(c(nodes[rand2[1],2], nodes[rand2[length(rand2)],2]) , c(nodes[rand2[1],3], nodes[rand2[length(rand2)],3]), type = "b", col = "blue", lty="dashed")

mylabels <- nodes[,1]
mylabels <- mylabels - 1
text(nodes[,2], nodes[,3], labels=mylabels, cex= 0.7, pos=3)

## PX cut 
# sg1 <- c(20, 47)
# sg2 <- c(13, 24)
# 
# lines(c(nodes[sg1[1]+1,2], nodes[sg1[2]+1,2]), c(nodes[sg1[1]+1,3], nodes[sg1[2]+1,3]),lty = "dotted", col = "green", lwd=3)
# lines(c(nodes[sg2[1]+1,2], nodes[sg2[2]+1,2]), c(nodes[sg2[1]+1,3], nodes[sg2[2]+1,3]),lty = "dotted", col = "green", lwd=3)


# 
# 
# # my solutions
# afterpx1 <- c(115, 105, 106, 107, 108, 109, 75, 74, 73, 76, 77, 72, 71, 78, 79, 80, 81, 87, 88, 89, 90, 104, 91, 92, 93, 94, 95, 96, 97, 98, 100, 102, 101, 103, 99, 63, 62, 52, 61, 60, 59, 56, 55, 54, 53, 51, 48, 50, 49, 47, 46, 45, 57, 58, 37, 38, 44, 43, 42, 41, 40, 39, 36, 35, 34, 67, 64, 65, 66, 68, 69, 27, 26, 28, 29, 30, 31, 32, 33, 25, 24, 23, 19, 22, 21, 20, 18, 17, 16, 15, 14, 13, 12, 11, 10, 6, 9, 7, 8, 0, 1, 2, 3, 4, 5, 84, 85, 83, 86, 82, 70, 121, 122, 120, 119, 123, 124, 125, 126, 110, 111, 112, 113, 114, 118, 116, 117)
# afterpx2 <- c(5, 13, 14, 18, 19, 22, 23, 24, 25, 26, 27, 69, 68, 67, 66, 65, 64, 59, 56, 55, 60, 54, 53, 61, 62, 103, 104, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 101, 102, 100, 99, 63, 52, 51, 48, 50, 49, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 57, 58, 36, 35, 34, 28, 29, 30, 31, 32, 33, 21, 20, 17, 15, 16, 12, 11, 10, 6, 9, 7, 8, 1, 0, 123, 124, 125, 126, 110, 111, 113, 114, 118, 117, 116, 115, 112, 106, 74, 105, 107, 108, 109, 75, 76, 73, 72, 77, 78, 71, 70, 120, 119, 122, 121, 79, 80, 81, 82, 86, 83, 85, 84, 2, 3, 4)
# afterpx1 <- afterpx1+1
# afterpx2 <- afterpx2+1
# 
# for(i in 2:length(afterpx1)){
#   points(c(nodes[afterpx1[i-1],2], nodes[afterpx1[i],2]) , c(nodes[afterpx1[i-1],3], nodes[afterpx1[i],3]), type = "b", col = "black", lty="solid", lwd=4)
# }
# points(c(nodes[afterpx1[1],2], nodes[afterpx1[length(afterpx1)],2]) , c(nodes[afterpx1[1],3], nodes[afterpx1[length(afterpx1)],3]), type = "b", col = "black", lty="solid", lwd=4)
# 
# 
# 
# for(i in 2:length(afterpx2)){
#   points(c(nodes[afterpx2[i-1],2], nodes[afterpx2[i],2]) , c(nodes[afterpx2[i-1],3], nodes[afterpx2[i],3]), type = "b", col = "purple", lty="solid", lwd=4)
# }
# points(c(nodes[afterpx2[1],2], nodes[afterpx2[length(afterpx2)],2]) , c(nodes[afterpx2[1],3], nodes[afterpx2[length(afterpx2)],3]), type = "b", col = "purple", lty="dashed")

