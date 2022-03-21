#' ModelInput class
#'
#' @name ModelInput-class
#' @rdname ModelInput-class
#' @exportClass ModelInput
#' @slot A.hat The Laplacian of a line graph built from the co-regulation graphs, where 
#' each node corresponds to a pair of analytes.
#' @slot node.wise.prediction graph of co-regulation results from IntLIM
#' @slot true.phenotypes data frame of true phenotypes given analyte pairs
#' and clinical covariates
#' @slot outcome.type "numeric" or "categorical"
#' @slot coregulation.graph the original co-regulation graph for the input data.
#' @slot line.graph the line graph of the input data.
methods::setClass(
  Class="ModelInput",
  representation(A.hat="matrix",
                 node.wise.prediction="matrix",
                 true.phenotypes="numeric",
                 coregulation.graph="matrix",
                 line.graph="matrix",
                 outcome.type="character",
                 modified.outliers="list")
)
#' PoolingFilter class
#'
#' @name PoolingFilter-class
#' @rdname PoolingFilter-class
#' @exportClass PoolingFilter
#' @slot filter The filter that maps each dimension of the input to one of k clusters.
#' @slot filter.type One of "mean", "min", "median", or "max".
#' @slot cluster.sizes A vector of cluster sizes (by number of dimensions mapping
#' @slot individual.filters A list of filters corresponding to each individual sample.
#' to the cluster)
methods::setClass(
  Class="PoolingFilter",
  representation(filter="matrix",
                 filter.type="character",
                 cluster.sizes="matrix",
                 individual.filters="list")
)
#' ModelResults class
#'
#' @name ModelResults-class
#' @rdname ModelResults-class
#' @exportClass ModelResults
#' @slot model.input An object of the modelInput class.
#' @slot pooling.filter A matrix that pools convolution results.
#' @slot iteration.tracking A data frame to track the iteration, weight, and error
#' values for each iteration of training.
#' @slot max.iterations Maximum number of iterations.
#' @slot convergence.cutoff Cutoff for convergence.
#' @slot learning.rate Learning rate used during training.
#' @slot activation.type Character value. Must be "sigmoid", "tanh", or "softmax".
#' @slot current.importance.weights Importance weights used in the current iteration.
#' @slot previous.weights Weights used in the previous iteration.
#' @slot current.gradient Gradient calculated for this iteration.
#' @slot weights.after.pooling Whether to include the weights after the pooling
#' layer (as opposed to before). Must be a boolean.
#' @slot outcome.prediction The prediction of the outcome.
#' @slot optimization.type One of "BGD", "momentum", "adagrad", or "adam"
#' @slot current.iteration Iteration (changes at each time step)
#' @slot previous.momentum Momentum value used in momentum optimization
#' @slot previous.update.vector Previous value used to update weights, used in ADAM
#' optimization
#' @slot sum.square.gradients Sum of squared gradients over iterations, used in
#' Adagrad optimization
methods::setClass(
  Class="ModelResults",
  representation(model.input="ModelInput",
                 pooling.filter="PoolingFilter",
                 iteration.tracking="data.frame",
                 max.iterations="numeric",
                 convergence.cutoff="numeric",
                 learning.rate="numeric",
                 activation.type="character",
                 previous.importance.weights="matrix",
                 current.importance.weights="matrix",
                 current.gradient="matrix",
                 previous.update.vector="matrix",
                 current.iteration="numeric",
                 outcome.prediction="numeric",
                 optimization.type="character",
                 sum.square.gradients="matrix",
                 previous.momentum="matrix")
)


