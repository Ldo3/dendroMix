#' Plot dendrogram of mixing measures
#'
#' @param dmm Output from \code{dendrogram_mixing}.
#' @param dim Dimension to plot.
#' @param atom_size Point size.
#' @param lwd Line width.
#'
#' @return A base R plot.
#' @export
plot_dendrogram_mixing <- function(dmm, dim=1, atom_size=1, lwd=1) {
     Gs = dmm$Gs
     heights = dmm$hc$height
     merge_pair = dmm$merge_pair
     K_bar = length(dmm$Gs)

     heights_cumsum <- c(0, heights)[-(K_bar + 1)]

     plot(NULL, xlim = c(0, max(heights_cumsum)),
          ylim = range(sapply(as.matrix(Gs[[1]]$thetas)[,dim], function(x) range(x))),
          xlab = "Heights", ylab = paste0("theta[", dim, "]"), main = "Dendrogram of mixing measures")

     for (step in seq_len(K_bar - 1)) {
          Gs[[step]]$thetas = as.matrix(Gs[[step]]$thetas)
          Gs[[step+1]]$thetas = as.matrix(Gs[[step+1]]$thetas)

          i_merge = merge_pair[step, 1]
          j_merge = merge_pair[step, 2]
          for (i in seq_len(j_merge - 1)) {
               y1 <- Gs[[step]]$thetas[i, dim]
               y2 <- Gs[[step+1]]$thetas[i, dim]
               x1 <- heights_cumsum[step]
               x2 <- heights_cumsum[step + 1]
               segments(x1, y1, x2, y1, col = adjustcolor("black", alpha.f = .8), lwd=lwd)
               segments(x2, y1, x2, y2, col = adjustcolor("black", alpha.f = .8), lwd=lwd)
          }

          y1 <- Gs[[step]]$thetas[j_merge, dim]
          y2 <- Gs[[step+1]]$thetas[i_merge, dim]
          x1 <- heights_cumsum[step]
          x2 <- heights_cumsum[step + 1]
          segments(x1, y1, x2, y1, col = adjustcolor("black", alpha.f = .8), lwd=lwd)
          segments(x2, y1, x2, y2, col = adjustcolor("black", alpha.f = .8), lwd=lwd)

          for (i in c((j_merge+1):(K_bar-step+1))) {
               if (step < K_bar - 1){
                    y1 <- Gs[[step]]$thetas[i, dim]
                    y2 <- Gs[[step+1]]$thetas[i-1, dim]
                    x1 <- heights_cumsum[step]
                    x2 <- heights_cumsum[step + 1]
                    segments(x1, y1, x2, y1, col = adjustcolor("black", alpha.f = .8), lwd=lwd)
                    segments(x2, y1, x2, y2, col = adjustcolor("black", alpha.f = .8), lwd=lwd)
               }
          }
     }
     for (i in seq_len(K_bar)) {
          points(rep(heights_cumsum[i], K_bar - i + 1), Gs[[i]]$thetas[, dim], pch = 19, cex = atom_size, col=i)
     }
     text(rep(heights_cumsum[1], K_bar), Gs[[1]]$thetas[, dim],
          labels = seq_len(K_bar), col = 'white', cex = atom_size / 3)
}
