library(keras)
# Even Neural Networks still can't simulate better
# 2018-07-31
predictors <- c("mean", "sd", "cv", "skewness", "kurtosis", "lat")
response   <- c("lambda")
vars       <- c(predictors, response)
npredict   <- length(predictors)


INPUT_lst <- znorm(temp[, ..vars])

data <- INPUT_lst$data
# data <- as.matrix(temp[, vars, with = F])
# data <- alply(data, 2, zscore) %>% do.call(cbind, .) %>% set_colnames(vars)
n    <- nrow(data)
I    <- sample(1:n, n*0.3)

INPUT <- list(
    train = data[-I, ],
    test  = data[ I, ])

x_train <- INPUT$train[, predictors]
y_train <- INPUT$train[, response, drop = F]
x_test  <- INPUT$test[, predictors]
y_test  <- INPUT$test[, response, drop = F]

# ------------------------------------------------------------------------
# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), npredict))
x_test  <- array_reshape(x_test, c(nrow(x_test), npredict))

# rescale
# x_train <- x_train / 255
# x_test <- x_test / 255

# ------------------------------------------------------------------------
# y_train <- to_categorical(y_train, 10)
# y_test  <- to_categorical(y_test, 10)

# ------------------------------------------------------------------------
model <- keras_model_sequential()
activation <- "linear"

model %>%
    layer_dense(units = 4, activation = activation, input_shape = c(npredict)) %>%
    # layer_dense(units = 4, activation = activation) %>%
    layer_dense(units = 1, activation = activation)
    # layer_dense(units = 4, activation = activation) # %>%
    # layer_dropout(rate = 0.3) %>%
    # layer_dropout(rate = 0.4) %>%
    # layer_dense(units = 10, activation = 'softmax')

# ------------------------------------------------------------------------
summary(model)

# ------------------------------------------------------------------------
model %>% compile(
    # loss = 'mse',
    # optimizer = optimizer_rmsprop()
    loss='mse',
    optimizer=optimizer_sgd(lr=0.01, clipnorm=0.5)
    # metrics = c('mse')
)

# ---- results='hide'-----------------------------------------------------
history <- model %>% fit(
    x_train, y_train,
    epochs = 30, batch_size = 2000,
    validation_split = 0.2
)

# ------------------------------------------------------------------------
plot(history)

# ---- results = 'hide'---------------------------------------------------
model %>% evaluate(x_test, y_test)

# ---- results = 'hide'---------------------------------------------------
z <- model %>% predict(x_test)

GOF(y_test, z)

# ggplot(data.frame(obs = y_test[, 1], sim = z), aes(obs, sim)) +
#     geom_point() +
#     geom_abline(slope = 1, color = "red") +
#     geom_density2d()
# plot(y_test, z); grid(); abline(a = 0, b = 1, col = "red")
