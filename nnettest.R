library(nnet)
n = 3
redwine = read.csv("winequality-red.csv", sep = ';')
nrowwine = nrow(redwine)
ncolwine = ncol(redwine)
redwine$quality = as.factor(redwine$quality)
if (n == 1)
  set.seed(123456)

tp_eachcorrects = NULL;
tp_allcorrects = NULL;
tp = seq(0.5,0.9,0.1)
for (tr in tp) {
  eachcorrects = numeric(length(levels(redwine$quality)))
  names(eachcorrects) = levels(redwine$quality)
  eachn = numeric(length(eachcorrects))
  names(eachn) = levels(redwine$quality)
  allcorrects = numeric(n)
  for (i in 1:n) {
    id = sample(1:nrowwine, nrowwine * tr)
    winetrain = redwine[id,]
    winetest = redwine[setdiff(1:nrowwine, id),]
    ideal = class.ind(winetrain$quality)
    winenet = nnet(winetrain[, -ncolwine], ideal, size = 30, softmax = TRUE, maxit = 20)
    winepredict = predict(winenet, winetest[,-ncolwine], 'class')
    cmat = table(winetest[,ncolwine], winepredict)
    
    eachcorrectrate = lapply(rownames(cmat), function(x) {
      if (x %in% colnames(cmat)) {
        cmat[x,x] / sum(cmat[x,])
      } else
        0
    })
    
    eachcorrects[rownames(cmat)] = eachcorrects[rownames(cmat)] + unlist(eachcorrectrate)
    eachn[rownames(cmat)] = eachn[rownames(cmat)] + 1
    
    allcorrectrate = lapply(rownames(cmat), function(x) {
      if (x %in% colnames(cmat)) {
        cmat[x,x]
      } else
        0
    })
    allcorrects[i] = sum(unlist(allcorrectrate)) / sum(cmat)
  }
  eachcorrects = eachcorrects / eachn
  tp_eachcorrects = rbind(tp_eachcorrects, eachcorrects)
  tp_allcorrects = rbind(tp_allcorrects, allcorrects)
}
rownames(tp_eachcorrects) = tp
rownames(tp_allcorrects) = tp
par(mfrow = c(2,2))

plot(rownames(tp_allcorrects), apply(tp_allcorrects, 1, mean), type = 'b', xlab = '训练比例', ylab = '平均总正确比例')
plot(apply(tp_allcorrects, 2, mean), type = 'b', xlab = '迭代次数', ylab = '平均总正确比例')

plot(colnames(tp_eachcorrects), apply(tp_eachcorrects, 2, mean, na.rm = TRUE), type = 'b', xlab = '类别', ylab = '不同训练比列的平均正确率')
