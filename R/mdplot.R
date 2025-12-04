mdplot <- function(delta, sd = 1, constant = 3, n = 501, alpha = 0.25, lwd = 1, 
                   col1 = "darkblue", col2 = "darkred", xlim, addSens = TRUE, 
                   digits = 3){
  stopifnot(!missing(delta))
  stopifnot(is.numeric(delta))
  stopifnot(is.numeric(sd))
  stopifnot(is.numeric(constant))
  n <- as.integer(n)
  if(length(delta) != 1 || length(sd) != 1 || 
     length(constant) != 1 || length(n) != 1){
    stop("'delta', 'sd', 'constant' and 'n' must be single numbers (vectors of length 1)!")
  }
  if(delta < 0){
    stop("'delta' must be a non-negative real number!")
  } 
  if(sd <= 0){
    stop("'sd' must be a positive real number!")
  } 
  
  if(addSens){
    sens <- pnorm(delta/2, mean = 0, sd = sd)
  }
  Min <- -constant*sd
  Max <- delta + constant*sd
  if(missing(xlim)) xlim <- c(Min, Max)
  DF <- data.frame(x = seq(from = Min, to = Max, length = n))
  gg <- ggplot(data = DF, aes(x = x)) + 
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
  print(gg)
}
