volcano <- function(x, ...){ 
  UseMethod("volcano")
}
volcano.default <- function(x, pval, effect0 = 0, sig.level = 0.05, 
                            effect.low = NULL, effect.high = NULL,
                            color.low = "#4575B4", color.high = "#D73027",
                            xlab = "effect size", ylab = "-log10(p value)",
                            title = "Volcano Plot", alpha = 1, shape = 19, ...){
  effect <- x
  stopifnot(is.numeric(effect), is.numeric(pval))
  stopifnot(all(pval > 0), all(pval <= 1))
  stopifnot(length(effect0) == 1)
  stopifnot(length(sig.level) == 1)
  stopifnot(sig.level > 0, sig.level < 1)
  if(!is.null(effect.low) & !is.null(effect.high))
    stopifnot(length(effect.low) == 1, length(effect.high) == 1, 
              effect.low < effect.high)
  stopifnot(length(color.low) == 1, length(color.high) == 1)
  
  if(!is.null(effect.low)){
    low <- effect < effect.low & pval < sig.level
  }else{
    low <- logical(length(effect))
  }
  if(!is.null(effect.high)){
    high <- effect > effect.high & pval < sig.level
  }else{
    high <- logical(length(effect))
  }
  effect.color <- ifelse(low, "low", "normal")
  effect.color[high] <- "high"
  effect.color <- factor(effect.color, levels = c("low", "normal", "high"))
  effect.color <- factor(effect.color)
  
  DF <- data.frame(effect = effect, pval = pval, col = effect.color)
  if(is.null(effect.low) & is.null(effect.high)){
    gg <- ggplot(DF, aes(x = effect, y = pval)) + 
      geom_point(alpha = alpha, shape = shape) + scale_y_neglog10() + 
      geom_vline(xintercept = effect0) + 
      geom_vline(xintercept = c(effect.low, effect.high), 
                 color = c(color.low, color.high)) +
      geom_hline(yintercept = sig.level) +
      labs(x = xlab, y = ylab) +
      ggtitle(title)
  }else{
    gg <- ggplot(DF, aes(x = effect, y = pval, color = col)) + 
      geom_point(alpha = alpha, shape = shape) + scale_y_neglog10() + 
      scale_color_manual(values = c("low" = color.low, 
                                    "normal" = "black", 
                                    "high" = color.high)) +
      geom_vline(xintercept = effect0) + 
      geom_vline(xintercept = c(effect.low, effect.high), 
                 color = c(color.low, color.high)) +
      geom_hline(yintercept = sig.level) +
      labs(x = xlab, y = ylab) +
      ggtitle(title)
  }
  gg
}
