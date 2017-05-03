library(nnet)
n = 3
redwine = read.csv("winequality-red.csv", sep = ';')
nrowwine = nrow(redwine)
ncolwine = ncol(redwine)
redwine$quality = as.factor(redwine$quality)
if (n == 1)
  set.seed(123456)

eachcorrects = numeric(length(levels(redwine$quality)))
names(eachcorrects) = levels(redwine$quality)
eachn = numeric(length(eachcorrects))
names(eachn) = levels(redwine$quality)
allcorrects = numeric(n)
for (i in 1:n) {
  id = sample(1:nrowwine, nrowwine * 0.7)
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
