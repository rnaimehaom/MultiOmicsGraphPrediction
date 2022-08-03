#' ModelInput class
#'
#' @name ModelInput-class
#' @rdname ModelInput-class
#' @exportClass ModelInput
#' @slot edge.wise.prediction graph of co-regulation results from IntLIM
#' @slot true.phenotypes data frame of true phenotypes given analyte pairs
#' and clinical covariates
#' @slot coregulation.graph the original co-regulation graph for the input data.
#' @slot metaFeatures A list of calculated meta-feature data frames for each sample
#' @slot model.properties A data frame that includes model information, i.e. R^2,
#' interaction term p-value, and coefficients.
#' @slot input.data An IntLIMData object that includes slots for the sample data,
#' the analyte data, and the analyte meta data.
#' @slot covariates A list of covariates to use in the models.
#' @slot stype Outcome/phenotype column name.
#' @slot stype.class Class of outcome/phenotype (either "numeric" or "character").
methods::setClass(
  Class="ModelInput",
  representation(edge.wise.prediction="matrix",
                 true.phenotypes="numeric",
                 coregulation.graph="matrix",
                 metaFeatures="list",
                 model.properties = "data.frame",
                 input.data = "IntLimData",
                 covariates = "list",
                 stype = "character",
                 stype.class = "character",
                 outcome = "numeric",
                 independent.var.type = "numeric")
)

#' ModelResults class
#'
#' @name ModelResults-class
#' @rdname ModelResults-class
#' @exportClass ModelResults
#' @slot model.input An object of the modelInput class.
#' @slot iteration.tracking A data frame to track the iteration, weight, and error
#' values for each iteration of training.
#' @slot max.iterations Maximum number of iterations.
#' @slot convergence.cutoff Cutoff for convergence.
#' @slot learning.rate Learning rate used during training.
#' @slot activation.type Character value. Must be "sigmoid", "tanh", or "softmax".
#' @slot current.metaFeature.weights Metafeature weights used in the current iteration.
#' @slot previous.weights Weights used in the previous iteration.
#' @slot current.gradient Gradient calculated for this iteration.
#' @slot optimization.type One of "SGD", "momentum", "adagrad", or "adam"
#' @slot current.iteration Iteration (changes at each time step)
#' @slot previous.momentum Momentum value used in momentum optimization
#' @slot previous.update.vector Previous value used to update weights, used in ADAM
#' optimization
#' @slot sum.square.gradients Sum of squared gradients over iterations, used in
#' Adagrad optimization
methods::setClass(
  Class="ModelResults",
  representation(model.input="ModelInput",
                 iteration.tracking="data.frame",
                 max.iterations="numeric",
                 convergence.cutoff="numeric",
                 learning.rate="numeric",
                 activation.type="character",
                 previous.metaFeature.weights="matrix",
                 current.metaFeature.weights="matrix",
                 current.gradient="matrix",
                 previous.update.vector="matrix",
                 current.iteration="numeric",
                 optimization.type="character",
                 sum.square.gradients="matrix",
                 previous.momentum="matrix",
                 pairs="character")
)


