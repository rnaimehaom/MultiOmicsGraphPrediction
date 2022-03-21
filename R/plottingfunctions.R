#' For each sample, save a plot of phenotype predictions from each edge in the
#' IntLIM graph. The edge weight is the phenotype prediction
#' using the pair of analytes connected to the edge. Weights correspond to color
#' in the graph. A color bar is shown for reference, and the true phenotype is
#' listed at the top of the plot.
#' @param graphWithPredictions An igraph object. This graph is the co-regulation graph
#' generated using IntLIM analysis of analyte pairs. Weights correspond to phenotype
#' predictions.
#' @param inputData Named list (output of FilterData()) with gene expression,
#' metabolite abundances, and associated meta-data
#' @param stype The phenotype of interest. This should correspond to a column in the
#' input data.
#' @param dirName The name of the directory where the output images will be saved.
#' @param continuous A boolean indicating whether the phenotype is continuous (TRUE)
#' or discrete (FALSE).
#' @export
SaveGraphPredictionPlots <- function(graphWithPredictions, inputData, stype, dirName,
                                     continuous = TRUE){
  # Create folder.
  dir.create(dirName)
  
  # Save each graph.
  for(name in names(graphWithPredictions)){
    
    # Connect to the file.
    grDevices::png(paste0(paste(dirName, make.names(name), sep = "\\"), ".png"))
    
    # Extract graph for the subject of interest.
    g <- graphWithPredictions[[name]]
    
    # Set up the layout and margins.
    graphics::layout(t(1:2),widths=c(5.5,1.5))
    graphics::par(mar=c(0,0,2,3))
    
    # Plot the graph.
    plot(g, layout = igraph::layout.fruchterman.reingold, vertex.label = NA)
    
    # Add the true phenotype and subject ID.
    inputDat <- inputData@phenoData$expression$main@data
    true_phen <- inputDat[which(rownames(inputDat) == name), stype]
    graphics::text(x = -1, y = 1.2, paste0(name, " (true phenotype is ", true_phen, ")"), pos = 4)
    
    # Add information about conversion to factors for discrete data.
    if(continuous == FALSE){
      inputDataPhen <- as.factor(inputDat[,stype])
      graphics::text(x = -1, y = 1.1, paste(levels(inputDataPhen)[1], " <= 0"), pos = 4)
      graphics::text(x = -1, y = 1.0, paste(levels(inputDataPhen)[2], " >= 1"), pos = 4)
    }
    
    # Add the color bar.
    color <- igraph::edge_attr(g, name = "color")[order(igraph::edge_attr(g, name = "weight"))]
    labs <- igraph::edge_attr(g, name = "weight")[order(igraph::edge_attr(g, name = "weight"))]
    lab_quants <- seq(min(labs), max(labs), by = (max(labs)-min(labs))/5)
    graphics::image(y=1:100,z=t(1:100), col=color, axes=FALSE, main="Prediction", cex.main=.8)
    graphics::axis(4,cex.axis=0.8, at = seq(0, 100, by = 20), labels = format(as.list(lab_quants), 
                                                                    digits=0, 
                                                                    scientific=FALSE), 
         las = 1)
    
    # Close the file connection.
    grDevices::dev.off()
  }
}

#' Plot the graph with edges colored by weight in the final outcome.
#' @param graph The co-regulation graph.
#' @param results A modelResults object.
#' @export
PlotGraphWeights <- function(graph, results){
  
  # Set up the layout and margins.
  graphics::layout(t(1:2),widths=c(5.5,1.5))
  graphics::par(mar=c(0,0,2,3))
  
  # Match the weights to graph edges.
  g <- igraph::as_data_frame(graph)
  g$to <- make.names(g$to)
  g$from <- make.names(g$from)
  weights <- results@current.weights
  if(results@weights.after.pooling == TRUE){
    S <- results@pooling.filter@filter
    weights <- t(matrix(rep(weights, dim(S)[1]), ncol = dim(S)[1]))
    sum_S <- colSums(S)
    S.weighted <- S * weights / sum_S
    S.flat <- rowSums(S.weighted)
    weights <- S.flat
  }
  names(weights) <- rownames(results@model.input@node.wise.prediction)
  weights_by_edge_name <- lapply(1:dim(g)[1], function(edge){
    forwards <- paste(g$from[edge], g$to[edge], sep = "__")
    backwards <- paste(g$to[edge], g$from[edge], sep = "__")
    which_weight <- union(which(names(weights) == forwards), 
                          which(names(weights) == backwards))
    the_weight <- NA
    if(length(which_weight) > 0){
      the_weight <- weights[which_weight]
    }
    return(the_weight)
  })
  
  # Add the weights to the data frame.
  g$weight <- unlist(weights_by_edge_name)
  g <- g[which(!is.na(g$weight)),]
  
  # Map weights to colors.
  pal <- grDevices::colorRampPalette(c("blue", "red"))(100)
  range_weight <- range(g$weight)
  color_scale <- pal[findInterval(g$weight, seq(range_weight[1], range_weight[2], 
                                            length.out = length(pal)+1), all.inside = TRUE)]
  g$color <- color_scale

  # Plot the graph.
  new_graph <- igraph::graph_from_data_frame(g, directed = FALSE)
  plot(new_graph, layout = igraph::layout.fruchterman.reingold, vertex.label = NA,
       vertex.size = 3)
  
  # Add the color bar.
  color <- igraph::edge_attr(new_graph, name = "color")[order(igraph::edge_attr(new_graph, 
                                                                                name = "weight"))]
  labs <- igraph::edge_attr(new_graph, name = "weight")[order(igraph::edge_attr(new_graph, 
                                                                                name = "weight"))]
  lab_quants <- seq(min(labs), max(labs), by = (max(labs)-min(labs))/5)
  graphics::image(y=1:100,z=t(1:100), col=color, axes=FALSE, main="Weight", cex.main=.8)
  graphics::axis(4,cex.axis=0.8, at = seq(0, 100, by = 20), labels = format(as.list(lab_quants), 
                                                                            digits=0, 
                                                                            scientific=FALSE), 
                 las = 1)
}

#' Plot the graph as a heatmap with edges colored by interaction coefficient.
#' @param inputResults The IntLIM results (from RunIntLim())
#' @param inputData The input data (from ReadData())
#' @export
PlotGraphWeightsHeatmap <- function(inputResults, inputData){
  # Add the analytes to the data frame.
  edge_df = data.frame(Analyte.1 = inputResults[,1], 
                       Analyte.2 = inputResults[,2])
  
  # Add the weights and corresponding colors.
  edge_df$Interaction.Coeff = inputResults$interaction_coeff
  
  # Truncate analyte names.
  edge_df$Analyte.1 <- unlist(lapply(edge_df$Analyte.1, function(a){
    return(substr(a, 1, 25))
  }))
  edge_df$Analyte.2 <- unlist(lapply(edge_df$Analyte.2, function(a){
    return(substr(a, 1, 25))
  }))
  
  # Set variables to NULL to appease R CMD check (these lines
  # serve no purpose other than that).
  Analyte.1 <- NULL
  Analyte.2 <- NULL
  Interaction.Coeff <- NULL
  
  # Plot
  plt <- ggplot2::ggplot(edge_df, ggplot2::aes(x=Analyte.1, 
                                               y=Analyte.2, 
                                               fill=Interaction.Coeff)) +
    ggplot2::scale_fill_gradient(low = "red", high = "blue") + ggplot2::geom_tile() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  print(plt)
}

#' A wrapper for PlotGraphWeights that applies the function across all folds.
#' @param inputResults The IntLIM results (from RunIntLimAllFolds())
#' @param inputData The input data (from CreateCrossValFolds())
#' @export
PlotGraphWeightsHeatmapAllFolds <- function(inputResults, inputData){
  return(lapply(1:length(inputResults), function(i){
    PlotGraphWeightsHeatmap(inputResults = inputResults[[i]], inputData = inputData[[i]])
  }))
}

#' Plot the graph with positive associations colored blue and negative associations
#' colored red.
#' @param graph The co-regulation graph.
#' @param saveInFile Location where the file should be saved. If NULL, then the output
#' is plotted without being saved. Default is NULL.
#' @param vertices List of vectors of vertices to plot. This is used if one
#' wishes to focus on a subset of vertices. If NULL, then all vertices are plotted.
#' Default is NULL.
#' @param truncateTo Vertex names are truncated to the first "truncateTo" characters.
#' Default is 4. If NULL, names are not truncated.
#' @param title Title of plot
#' @export
PlotCoRegulationGraph <- function(graph, title, saveInFile = NULL, vertices = NULL,
                                  truncateTo = 4){
  
  # Extract subgraph.
  if(!is.null(vertices) && !is.null(vertices)){
    vert_id <- lapply(vertices, function(v){
      return(match(v, igraph::V(graph)$name))
    })
    graph <- igraph::subgraph(graph, v = vert_id)
  }
  
  # Truncate names.
  graph_labels <- igraph::V(graph)$name
  if(!is.null(truncateTo)){
    graph_labels <- unlist(lapply(igraph::V(graph)$name, function(v){
      return(paste0(substr(v, 1, truncateTo), "."))
    }))
  }
  
  # Save or plot.
  if(is.null(saveInFile)){
    plot(graph, main = title, layout = igraph::layout.random, 
         vertex.label = graph_labels)
  }else{
    grDevices::png(saveInFile,res = 1200)
    plot(graph, main = title, layout = igraph::layout.random, 
         vertex.label = graph_labels)
    grDevices::dev.off()
    print(paste("Saved plot to", saveInFile))
  }
}

#' Wrapper for PlotCoRegulationGraph.
#' @param graph The co-regulation graph.
#' @param saveInDir Directory where file should be saved. If NULL, then the output
#' is plotted without being saved. Default is NULL
#' @param vertices List of vectors of vertices to plot. This is used if one
#' wishes to focus on a subset of vertices. If NULL, then all vertices are plotted.
#' Default is NULL
#' @param truncateTo Vertex names are truncated to the first "truncateTo" characters.
#' Default is 4. If NULL, names are not truncated.
#' @export
PlotCoRegulationGraphAllFolds <- function(graph, saveInDir = NULL, vertices = NULL,
                                  truncateTo = 4){
  for(i in 1:length(graph)){
    g <- graph[[i]]
    if(!is.null(saveInDir)){
      fileName <- paste(saveInDir, paste0("coreg_graph", i, ".png"))
    }
    PlotCoRegulationGraph(graph=g, saveInFile=saveInDir, vertices=vertices,
                          truncateTo=truncateTo, title=paste("Fold", i))
  }
}

#' Plot the graph for each sample with edges colored according to prediction.
#' Include a color scale for the predictions.
#' @param graph The graph with predictions projected onto it.
#' @param inputData Named list (output of 
#' FilterData()) with gene expression, metabolite abundances, 
#' and associated meta-data
#' @param saveInDir Directory where files should be saved. If NULL, then the output
#' is plotted without being saved. Default is NULL.
#' @param vertices List of vectors of vertices to plot. This is used if one
#' wishes to focus on a subset of vertices. If NULL, then all vertices are plotted.
#' Default is NULL.
#' @param truncateTo Vertex names are truncated to the first "truncateTo" characters.
#' Default is 4. If NULL, names are not truncated.
#' @param includeLabels whether or not to include labels. Defaults to TRUE.
#' @param cutoffs Cutoff weight values, which can be included for visibility.
#' Default is NULL.
#' @param vertexSize Vertex size to use in display.
#' @export
PlotGraphPredictions <- function(graph, inputData, saveInDir = NULL, 
                                 vertices = NULL, truncateTo = 4, includeLabels = TRUE,
                                 cutoffs = NULL, vertexSize = 10){
  for(j in 1:length(graph)){
    g <- graph[[j]]
    
    # Extract subgraph.
    if(!is.null(vertices) && !is.null(vertices)){
      vert_id <- lapply(vertices, function(v){
        return(match(v, igraph::V(g)$name))
      })
      g <- igraph::subgraph(g, v = vert_id)
    }
    
    # Set up variables for plotting.
    bin_count <- 100
    title <- paste(names(graph)[j], paste("True Outcome =",
                                          formatC(inputData$p[[j]], digits = 2, 
                                                  format = "f")), sep = "\n")
    edges <- igraph::E(g)
    wt <- edges$weight
    wt[which(wt<(-1))]<-cutoffs[1]
    wt[which(wt>1)]<-cutoffs[2]
    intervals <- seq(range(wt)[1], range(wt)[2], by = (range(wt)[2] - range(wt)[1]) / (bin_count - 1))
    subject_color_scale <- findInterval(wt, intervals)
    pal <- grDevices::colorRampPalette(c("limegreen", "purple"))(bin_count+1)
    igraph::E(g)$color <- pal[subject_color_scale]
    minPred <- min(wt)
    maxPred <- max(wt)
    ticks <- seq(minPred, maxPred, len=11)
    scale <- (bin_count + 1)/(maxPred-minPred)
    colorBarTitle <- "Prediction"
    
    # Plot.
    graph_labels = NA
    if(includeLabels == TRUE){
      graph_labels <- unlist(lapply(igraph::V(g)$name, function(v){
        return(paste0(substr(v, 1, truncateTo), "."))
      }))
    }

    graphics::par(mfrow=c(1,2))
    plot(g, main = title,
         vertex.label = graph_labels, vertex.frame.color = "black", vertex.size = vertexSize,
         edge.arrow.size = 0.5)
    plot(c(0,10), c(minPred,maxPred), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
    title(colorBarTitle, adj = 0, cex.main = 0.75)
    graphics::axis(2, ticks, las=1)
    for (l in 1:(length(pal)-1)) {
      y <- (l-1)/scale + minPred
      graphics::rect(0,y,1,y+1/scale, col=pal[l], border=NA)
    }
  }
}

#' Wrapper for PlotGraphPredictions.
#' @param graphs The graph with predictions projected onto it.
#' @param inputDataFolds List of named lists (output of 
#' CreateCrossValidationFolds()) with gene expression, metabolite abundances, 
#' and associated meta-data
#' @param saveInDir Directory where file should be saved. If NULL, then the output
#' is plotted without being saved. Default is NULL.
#' @param vertices List of vectors of vertices to plot. This is used if one
#' wishes to focus on a subset of vertices. If NULL, then all vertices are plotted.
#' Default is NULL.
#' @param truncateTo Vertex names are truncated to the first "truncateTo" characters.
#' Default is 4. If NULL, names are not truncated.
#' @export
PlotGraphPredictionsAllFolds <- function(graphs, inputDataFolds,
                                         saveInDir = NULL, vertices = NULL,
                                         truncateTo = 4){
  for(i in 1:length(graphs)){
    # Create temporary directory.
    saveInDir_tmp <- saveInDir
    if(!is.null(saveInDir)){
      saveInDir_tmp <- paste(saveInDir, paste0("fold",i), sep = "//")
    }
    PlotGraphPredictions(graph=graphs[[i]], inputDataFolds[[i]]$training,
                         saveInDir=saveInDir_tmp,
                         vertices=vertices, truncateTo=truncateTo)
  }
}

#' Plot the line graph for each sample with nodes colored according to prediction.
#' Include a color scale for the predictions.
#' @param modelInput A list of ModelInput objects.
#' @param saveInDir Directory where file should be saved. If NULL, then the output
#' is plotted without being saved. Default is NULL.
#' @param stype Outcome / phenotype
#' @param analytes List of vectors of analytes to plot. This is used if one
#' wishes to focus on a subset of analytes If NULL, then all vertices are plotted.
#' Default is NULL.
#' @param truncateTo Analyte names are truncated to the first "truncateTo" characters.
#' Default is 4. If NULL, names are not truncated.
#' @param weights The weight assigned to each node. This is encoded using opacity.
#' If NULL, all nodes are opaque. Default is NULL.
#' @param analytes analytes of interest. Default is NULL
#' @param includeLabels whether or not to include labels. Defaults to TRUE.
#' @param cutoffs Cutoff weight values, which can be included for visibility.
#' Default is NULL.
#' @param vertexSize Vertex size to use in display.
#' @export
PlotLineGraph <- function(modelInput, stype, saveInDir = NULL,
                          truncateTo = 2, weights = NULL,
                          analytes = NULL, includeLabels = TRUE, cutoffs = NULL,
                          vertexSize = 10){
  # Extract graph and predictions.
  line.graph <- modelInput@line.graph
  node.wise.prediction <- modelInput@node.wise.prediction
  Y <- modelInput@true.phenotypes
  if(length(unique(Y)) == 2 && min(Y) == 1){
    Y <- Y - 1
  }
  
  for(j in 1:dim(node.wise.prediction)[2]){
    
    # Build node and edge graphs.
    edge_df <- reshape2::melt(line.graph)
    edge_df <- edge_df[which(edge_df[,3] != 0),]
    edge_df$arrow.size <- 0.25
    color <- "black"
    node_df <- data.frame(names(node.wise.prediction[,j]), node.wise.prediction[,j],
                          color)
    rownames(node_df) <- names(node.wise.prediction[,j])
    colnames(node_df) <- c("name", "prediction", "color")
    node_df$frame.color <- color
    final_graph <- igraph::graph_from_data_frame(edge_df, vertices = node_df)
    
    # Filter graph if required.
    wt <- node.wise.prediction[,j]
    if(!is.null(analytes)){
      
      vertices_with_analytes <- unlist(lapply(analytes, function(a){
        return(igraph::V(final_graph)$name[grepl(a, igraph::V(final_graph)$name,
                                                 fixed = TRUE)])
      }))
      final_graph <- igraph::induced_subgraph(final_graph, vertices_with_analytes)
      wt <- node.wise.prediction[names(igraph::V(final_graph)),j]
    }
    wt[which(wt<cutoffs[1])]<-cutoffs[1]
    wt[which(wt>cutoffs[2])]<-cutoffs[2]
    
    # Set up node colors.
    bin_count <- 100
    # Make sure the spacing is even. We need to do this using seq.
    intervals <- seq(range(wt)[1], range(wt)[2],
                     by = (range(wt)[2] - range(wt)[1]) / (bin_count - 1))
    subject_color_scale <- findInterval(wt, intervals)
    pal <- grDevices::colorRampPalette(c("limegreen", "purple"))(bin_count+1)
    color <-pal[subject_color_scale]
    
    # Adjust color opacity.
    if(!is.null(weights)){
      opacity <- abs(weights) / max(abs(weights))
      color <- unlist(lapply(1:length(color), function(c){
        return(grDevices::adjustcolor(color[c], alpha.f = opacity[c]))
      }))
    }
    igraph::V(final_graph)$color <- color
    
    # Set up variables for plotting.
    title <- paste(colnames(node.wise.prediction)[j],
                   paste("True Outcome =", formatC(Y[j], digits = 2, format = "f")),
                   sep = "\n")
    minPred <- min(unname(wt))
    maxPred <- max(unname(wt))
    ticks <- seq(minPred, maxPred, len=11)
    scale <- (length(pal)-1)/(maxPred-minPred)
    colorBarTitle <- "Prediction"
    
    # Plot.
    graph_labels <- NA
    if(includeLabels == TRUE){
      graph_labels <- unlist(lapply(igraph::V(final_graph)$name, function(v){
        pieces <- strsplit(v, "__")[[1]]
        sub_from <- pieces[1]
        sub_to <- pieces[2]
        if(!is.null(truncateTo)){
          sub_from <- substr(pieces[1], 1, truncateTo)
          sub_to <- substr(pieces[2], 1, truncateTo)
        }
        return(paste(sub_from, sub_to, sep = "_"))
      }))
    }
    
    graphics::par(mfrow=c(1,2))
    plot(final_graph, main = title, vertex.label = graph_labels, vertex.size = vertexSize)
    plot(c(0,10), c(minPred,maxPred), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
    title(colorBarTitle, adj = 0, cex.main = 0.75)
    graphics::axis(2, ticks, las=1)
    for (l in 1:(length(pal)-1)) {
      y <- (l-1)/scale + minPred
      graphics::rect(0,y,1,y+1/scale, col=pal[l], border=NA)
    }
  }
}

#' Wrapper for PlotLineGraph.
#' @param modelInputs A list of ModelInput objects.
#' @param saveInDir Directory where file should be saved. If NULL, then the output
#' is plotted without being saved. Default is NULL.
#' @param stype Outcome / phenotype
#' @param analytes List of vectors of analytes to plot. This is used if one
#' wishes to focus on a subset of analytes If NULL, then all vertices are plotted.
#' Default is NULL.
#' @param truncateTo Analyte names are truncated to the first "truncateTo" characters.
#' Default is 2. If NULL, names are not truncated.
#' @param weights The list of weights assigned to each node. This is encoded using opacity.
#' If NULL, all nodes are opaque. Default is NULL.
#' @export
PlotLineGraphAllFolds <- function(modelInputs, stype, saveInDir = NULL, analytes = NULL,
                                         truncateTo = 2, weights=NULL){
  # Plot the line graphs.
  for(i in 1:length(modelInputs)){
    
    # Set up output directory.
    saveInDir_tmp <- saveInDir
    if(!is.null(saveInDir)){
      saveInDir_tmp <- paste(saveInDir, paste0("fold",i), sep = "//")
    }
    
    # Set up weights.
    wt <- weights
    if(!is.null(wt)){
      wt <- weights[[i]]
    }
    
    # Plot.
    PlotLineGraph(modelInput=modelInputs[[i]], stype=stype, saveInDir=saveInDir_tmp,
                  analytes=analytes[[i]], truncateTo=truncateTo, weights=wt)
  }
}

#' Plots a color-coded dendrogram of predictions for a given sample. The
#' dendrogram is defined by the line graph, and the colors correspond to the
#' prediction made by each pair.
#' @param modelInput A ModelInput object.
#' @param hierarchicalClustering The output of the doHierarchicalClustering
#' function.
#' @param sampleIndex Index of sample.
#' @param predictionLimits A vector of the upper and lower limits to plot.
#' Predictions falling outside of this range will be modified to fall within
#' this range. This is only for the purposes of plotting, so that outlier
#' predictions do not obscure differences within the rest of the plot.
#' @export
PlotPredictionDendrogram <- function(modelInput, hierarchicalClustering, sampleIndex, 
                                     predictionLimits){
  utils::globalVariables(c("yend", "label", "x", "y", "xend", "y.x", "pred"))
  
  # Process dendrogram.
  dendrogram_ends <- list()
  dendrogram_segments <- list()
  dend <- stats::as.dendrogram(hierarchicalClustering)
  dendrogram_data <- ggdendro::dendro_data(dend)
  dendrogram_segments <- dendrogram_data$segments
  dendrogram_ends <- dplyr::filter(dendrogram_segments, yend == 0)
  dendrogram_ends <- dplyr::left_join(dendrogram_ends, dendrogram_data$labels, by = "x")
  dendrogram_ends <- dplyr::rename(dendrogram_ends, node_name = label)
  
  # Plot prediction dendrogram.
  all_wt <- dendrogram_ends
  preds <- modelInput@node.wise.prediction[,sampleIndex]
  names(preds) <- rownames(modelInput@node.wise.prediction)
  all_wt$pred <- preds[dendrogram_ends$node_name]
  all_wt$pred[which(all_wt$pred > predictionLimits[2])] <- predictionLimits[2]
  all_wt$pred[which(all_wt$pred < predictionLimits[1])] <- predictionLimits[1]
  
  # Make the plot.
  ggplot2::ggplot() + 
    ggplot2::geom_segment(data = dendrogram_segments, 
                          ggplot2::aes_(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::geom_segment(data = all_wt, 
                          ggplot2::aes_(x=x, y=y.x, xend=xend, yend=yend, color=pred)) +
    ggplot2::scale_color_gradient(low = "#AAFF00", high = "#023020")+
    ggplot2::scale_y_reverse() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(), 
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(), 
                   axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(), 
                   axis.text.x=ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank(), 
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.ticks.y=ggplot2::element_blank()) +
    ggplot2::labs(color='Prediction') +
    ggplot2::ggtitle(paste(colnames(modelInput@node.wise.prediction)[sampleIndex], "-", 
                  "True Phenotype is",
                  formatC(modelInput@true.phenotypes[sampleIndex], digits = 2)))
}

#' Plots a dendrogram of the optimal subspace clustering.
#' @param optimalClustering The output of FindOptimalSubspaceClustering.
#' @export
PlotSubspaceClusteringDendrogram <- function(inputData,
                                             eigStep = 10, alphaMin = 0,
                                             alphaMax = 1, alphaStep = 0.1){
    
  # Find the optimal projection for best cluster separability.
  type1sim <- ComputeCosineSimilarity(t(inputData@analyteType1))
  type2sim <- ComputeCosineSimilarity(t(inputData@analyteType2))
  opt <- FindOptimalSubspaceClustering(type1Similarity = type1sim, 
                                       type2Similarity = type2sim,
                                       eigStep = eigStep, alphaMin = alphaMin,
                                       alphaMax = alphaMax, alphaStep = alphaStep)
  # Plot the dendrogram.
  plot(opt$dendrogram)
}

#' Plots a heatmap of the optimal subspace clustering.
#' @include 03_ComputeMetaFeatures.R
#' @param optimalClustering The output of FindOptimalSubspaceClustering.
#' @export
PlotSubspaceClusteringHeatmap <- function(inputData,
                                          eigStep = 10, alphaMin = 0,
                                          alphaMax = 1, alphaStep = 0.1){
  # Find the optimal projection for best cluster separability.
  type1sim <- ComputeCosineSimilarity(t(inputData@analyteType1))
  type2sim <- ComputeCosineSimilarity(t(inputData@analyteType2))
  opt <- FindOptimalSubspaceClustering(type1Similarity = type1sim, 
                                       type2Similarity = type2sim,
                                       eigStep = eigStep, alphaMin = alphaMin,
                                       alphaMax = alphaMax, alphaStep = alphaStep)
  
  # Plot the heatmap.
  stats::heatmap(ComputeCosineSimilarity(opt$L_mod), 
          Rowv = stats::as.dendrogram(opt$dendrogram), 
          Colv = stats::as.dendrogram(opt$dendrogram))
}
