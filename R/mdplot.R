mdplot <- function(delta, sd1 = 1, sd2 = 1, constant = 3, n = 501, alpha = 0.25, 
                    lwd = 1, col1 = "darkblue", col2 = "darkred", xlim, 
                    addSens = TRUE, digits = 3){
  stopifnot(!missing(delta))
  stopifnot(is.numeric(delta))
  stopifnot(is.numeric(sd1))
  stopifnot(is.numeric(sd2))
  stopifnot(is.numeric(constant))
  n <- as.integer(n)
  if(length(delta) != 1 || length(sd1) != 1 || length(sd2) != 1 ||
     length(constant) != 1 || length(n) != 1){
    stop("'delta', 'sd1', 'sd2', 'constant' and 'n' must be single numbers (vectors of length 1)!")
  }
  if(delta < 0){
    stop("'delta' must be a non-negative real number!")
  } 
  if(sd1 <= 0 || sd2 <= 0){
    stop("The values of 'sd1' and 'sd2' must be positive real numbers!")
  }
  
  var.equal <- ifelse(abs(sd1-sd2) < 1e-8, TRUE, FALSE)
  
  if(var.equal){
    sd <- sd1
    if(addSens){
      sens <- pnorm(delta/2, mean = 0, sd = sd)
    }
    Min <- -constant*sd
    Max <- delta + constant*sd
    if(missing(xlim)) xlim <- c(Min, Max)
    DF <- data.frame(x = seq(from = Min, to = Max, length = n))
    gg <- ggplot(data = DF, aes(x = .data$x)) + 
      stat_function(fun = dnorm, n = n, args = list(mean = 0, sd = sd), 
                    color = col1, lwd = lwd) +
      geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = sd),
                fill = col1, alpha = alpha) + 
      geom_vline(xintercept = 0, color = col1, lwd = lwd) +
      stat_function(fun = dnorm, n = n, args = list(mean = delta, sd = sd),
                    color = col2, lwd = lwd) +
      geom_area(stat = "function", fun = dnorm, args = list(mean = delta, sd = sd),
                fill = col2, alpha = alpha) +
      geom_vline(xintercept = delta, color = col2, lwd = lwd) +
      ylab("") + xlab("") + scale_y_continuous(breaks = NULL) +
      annotate(x = -Inf, y = Inf, geom = "text", 
               label = paste0("MD = ", delta, ", SD = ", sd),
               hjust = -0.15, vjust = 1.5, fontface = "bold") +
      ggtitle("Mean Difference (MD)")
    if(addSens){
      gg <- gg + geom_vline(xintercept = delta/2, lwd = lwd, linetype = "dashed") +
        annotate(x = Inf, y = Inf, geom = "text", 
                 label = paste0("sensitivity = specificity = ", signif(sens, digits)),
                 hjust = 1.15, vjust = 1.5, fontface = "bold") +
        scale_x_continuous(breaks = c(0, delta/2, delta), 
                           labels = c("0", expression(paste(frac(MD,2))), "MD"), limits = xlim)
    }else{
      gg <- gg + scale_x_continuous(breaks = c(0, delta), labels = c("0", "MD"), 
                                    limits = xlim)
    }
  }else{
    Min <- min(-constant*sd1, delta - constant*sd2)
    Max <- max(constant*sd1, delta + constant*sd2)
    if(addSens){
      diff <- function(x, delta, sd1, sd2){
        dnorm(x, mean = 0, sd = sd1) - dnorm(x, mean = delta, sd = sd2)
      }
      if(sd1 < sd2){
        cutoff <- uniroot(diff, interval = c(0, Max), 
                          delta = delta, sd1 = sd1, sd2 = sd2)$root
      }else{
        cutoff <- uniroot(diff, interval = c(Min, delta), 
                          delta = delta, sd1 = sd1, sd2 = sd2)$root
      }
      sens <- pnorm(cutoff, mean = 0, sd = sd1)
      spec <- pnorm(cutoff, mean = delta, sd = sd2, lower.tail = FALSE)
    }
    if(missing(xlim)) xlim <- c(Min, Max)
    DF <- data.frame(x = seq(from = Min, to = Max, length = n))
    gg <- ggplot(data = DF, aes(x = .data$x)) + 
      stat_function(fun = dnorm, n = n, args = list(mean = 0, sd = sd1), 
                    color = col1, lwd = lwd) +
      geom_area(stat = "function", fun = dnorm, args = list(mean = 0, sd = sd1),
                fill = col1, alpha = alpha) + 
      geom_vline(xintercept = 0, color = col1, lwd = lwd) +
      stat_function(fun = dnorm, n = n, args = list(mean = delta, sd = sd2),
                    color = col2, lwd = lwd) +
      geom_area(stat = "function", fun = dnorm, args = list(mean = delta, sd = sd2),
                fill = col2, alpha = alpha) +
      geom_vline(xintercept = delta, color = col2, lwd = lwd) +
      ylab("") + xlab("") + scale_y_continuous(breaks = NULL) +
      annotate(x = -Inf, y = Inf, geom = "text", 
               label = paste0("MD = ", delta, ", SD1 = ", sd1, ", SD2 = ", sd2),
               hjust = -0.15, vjust = 1.5, fontface = "bold") +
      ggtitle("Mean Difference (MD)")
    if(addSens){
      gg <- gg + geom_vline(xintercept = cutoff, lwd = lwd, linetype = "dashed") +
        annotate(x = Inf, y = Inf, geom = "text", 
                 label = paste0("sensitivity = ", signif(sens, digits), 
                                "\nspecificity = ", signif(spec, digits)),
                 hjust = 1.15, vjust = 1.5, fontface = "bold") +
        scale_x_continuous(breaks = c(0, cutoff, delta), 
                           labels = c("0", "cutoff", "MD"), limits = xlim)
    }else{
      gg <- gg + scale_x_continuous(breaks = c(0, delta), labels = c("0", "MD"), 
                                    limits = xlim)
    }
  }
  gg
}
