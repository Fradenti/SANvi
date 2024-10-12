# cl_col <- as.numeric(factor(apply(x$sim$RHO,1,which.max)))
# cl_row <- unique(unlist(lapply(x$sim$XI,
#                                function(x) apply(x, 1, which.max))))
# 
# boxplot(x$params$y~x$params$group,col=scales::alpha(cl_col, .5), xlab="Group",ylab="y",
#         main = paste(x$model, "- DC and posterior means"))
# abline(h = x$sim$theta_l[unique(cl_row),1],col=4,lty=2)
# 