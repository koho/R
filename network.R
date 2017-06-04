network = function(hiddenSizes, trainFcn = "traingd", ...) {
  numLayers = length(hiddenSizes) + 2
  extras = list(...)
  transferFcn = extras$transferFcn
  if (is.null(transferFcn))
    transferFcn = rep("sigmoid", numLayers - 1)
  else if (length(transferFcn) != numLayers - 1)
    stop("Transfer function length must agree")
  performFcn = if(is.null(extras$performFcn)) "mse" else extras$performFcn
  structure(list(numLayers = numLayers, sizes = c(NA, hiddenSizes, NA),
                 biases = NULL, weights = NULL,
                 trainFcn = trainFcn, performFcn = performFcn,
                 transferFcn = transferFcn),
            class = "network")
}

train = function(net, X, Y, ...) {
  UseMethod("train")
}

train.network = function(net, X, Y, epochs = 1000, goal = 0, ...) {
  net$sizes[c(1, net$numLayers)] = c(ncol(X), ncol(Y))
  net$biases = lapply(net$sizes[-1], function(x) matrix(rnorm(x)))
  net$weights = mapply(function(x, y) matrix(rnorm(x * y), x, y), 
                       net$sizes[-1], net$sizes[-net$numLayers], SIMPLIFY = FALSE)
  stopifnot(nrow(X) == nrow(Y))
  testx = list(...)$testx
  testy = list(...)$testy
  cat("\n", "Call:", "\n", sep = "")
  print(match.call(call = sys.call(-1)))
  cat("\nTraining Function: ", net$trainFcn, "\n", sep = "")
  cat("Training Set: ", paste(dim(X), collapse = "x"), "\n", sep = "")
  cat("Goal: ", as.character(goal), "\n\n", sep = "")
  progshow = epochs %/% (proglen <- 30)
  perfs = numeric()
  evaluation = numeric()
  for(epoch in seq_len(epochs)) {
    net = getFunction(net$trainFcn)(net, X, Y, ...)
    p = predict(net, X)
    perf = getFunction(net$performFcn)(net)$value(p, Y)
    perfs[epoch] = perf
    if(!is.null(testx) && !is.null(testy)) {
      test = predict(net, testx)
      evaluation[epoch] = sum(apply(test, 1, which.is.max) == apply(testy, 1, which.is.max)) / nrow(testx)
    }
    cat(paste0("\r", "[", paste0(rep("=", min(epoch %/% progshow, proglen)), collapse = ""), ">",
            paste0(rep(" ", proglen - min(epoch %/% progshow, proglen)), collapse = ""), "]",
            as.character(epoch), "/", as.character(epochs), "\t", "Perf: ",  as.character(perf)))
    if(perf <= goal) {
      cat("\n\nPerformance goal met.\n")
      break
    }
  }
  if(epoch >= epochs) cat("\n\nMaximum epoch reached.\n")
  net$train = list(epochs = epochs, goal = goal, perfs = perfs, test = evaluation, ...)
  net
}

traingd = function(net, X, Y, lr = 0.01,
                   batchSize = 10, ...) {
  n = nrow(X)
  trainInd = sample.int(n)
  miniBatches = lapply(seq(1, n, batchSize), function(x) x:min(x + batchSize - 1, n))
  for(b in seq_along(miniBatches)) {
    bInd = trainInd[miniBatches[[b]]]
    updates = updatebatch(net, X[bInd, ,drop = FALSE], Y[bInd, ,drop = FALSE], lr)
    net$biases = updates$biases
    net$weights = updates$weights
  }
  net
}

updatebatch = function(net, X, Y, lr) {
  nabla_b = lapply(lapply(net$biases, dim), function(size) array(0, size))
  nabla_w = lapply(lapply(net$weights, dim), function(size) array(0, size))
  n = nrow(X)
  X = t(as.matrix(X))
  Y = t(as.matrix(Y))
  for(i in seq_len(n)) {
    delta_nabla = backprop(net, X[,i,drop = FALSE], Y[,i,drop = FALSE])
    nabla_b = mapply(`+`, nabla_b, delta_nabla$nabla_b, SIMPLIFY = FALSE)
    nabla_w = mapply(`+`, nabla_w, delta_nabla$nabla_w, SIMPLIFY = FALSE)
  }
  weights = mapply(function(w, nw) w - (lr/n) * nw, net$weights, nabla_w, SIMPLIFY = FALSE)
  biases = mapply(function(b, nb) b - (lr/n) * nb, net$biases, nabla_b, SIMPLIFY = FALSE)
  list(weights = weights, biases = biases)
}

backprop = function(net, x, y) {
  nabla_b = lapply(lapply(net$biases, dim), function(size) array(0, size))
  nabla_w = lapply(lapply(net$weights, dim), function(size) array(0, size))
  # feedforward
  activation = x
  activations = vector("list", net$numLayers)
  activations[[1]] = x
  zs = vector("list", net$numLayers - 1)
  for(i in seq_along(net$weights)) {
    z = net$weights[[i]] %*% activation + net$biases[[i]]
    zs[[i]] = z
    activation = getFunction(net$transferFcn[i])(z)
    activations[[i + 1]] = activation
  }
  # backward pass
  #delta = (activation - y) * getFunction(paste0(net$transferFcn[net$numLayers - 1], "_prime"))(
  #  zs[[net$numLayers - 1]]
  #)
  delta = getFunction(net$performFcn)(net)$delta(zs[[net$numLayers - 1]], activation, y)
  nabla_b[[net$numLayers - 1]] = delta
  nabla_w[[net$numLayers - 1]] = delta %*% t(activations[[net$numLayers - 1]])
  for(l in seq(net$numLayers - 2, 1)) {
    delta = (t(net$weights[[l + 1]]) %*% delta) * getFunction(paste0(net$transferFcn[l], "_prime"))(
      zs[[l]]
    )
    nabla_b[[l]] = delta
    nabla_w[[l]] = delta %*% t(activations[[l]])
  }
  list(nabla_b = nabla_b, nabla_w = nabla_w)
}

predict.network = function(net, X) {
  n = nrow(X)
  X = t(as.matrix(X))
  Y = matrix(nrow = n, ncol = net$sizes[net$numLayers])
  for(i in seq_len(n)) {
    a = X[, i, drop = FALSE]
    for(l in seq_len(net$numLayers - 1)) {
      a = getFunction(net$transferFcn[l])(
        (net$weights[[l]] %*% a) + net$biases[[l]]
      )
    }
    Y[i,] = t(a)
  }
  Y
}

showPerformance = function(net, ...) {
  plot(net$train$perfs, type = "l", col = "blue", xlab = "Epoch", ylab = net$performFcn, main = "Performance", ...)
  abline(net$train$goal, 0, col = "red")
}

plot.network = function(net, ...) {
  showPerformance(net, ...)
}

print.network = function(net) {
  cat("\nLayers:", net$sizes, "\n")
  cat("Training Function:", net$trainFcn, "\n")
  cat("Cost Function:", net$performFcn, "\n")
  cat("Transfer Function:", net$transferFcn, "\n\n")
}

crossentropy = function(net) {
  list(
    value = function(a, y) sum(-y * log(a) - (1 - y) * log(1 - a), na.rm = TRUE) / nrow(a),
    delta = function(z, a, y) a - y
  )
}

mse = function(net) {
  list(
    value = function(a, y) sum(apply((a - y) ^ 2, 1, sum)) / nrow(a),
    delta = function(z, a, y) (a - y) * getFunction(paste0(net$transferFcn[net$numLayers - 1], "_prime"))(z)
  )
}

sigmoid = function(z) 1 / (1 + exp(-z))

sigmoid_prime = function(z) sigmoid(z) * (1 - sigmoid(z))

purelin = function(z) z

purelin_prime = function(z) 1

softmax = function(z) exp(z) / sum(exp(z))

softmax_prime = function(z) {}
