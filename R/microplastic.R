#' R package for microplastic image analysis
#' get_contour_plot
#'
#' @author Junah Bang \email{juna0033@@g.skku.edu} and Youngdeok Hwang \email{yhwang@@g.skku.edu}
#' @param img image source
#' @param min_val threshold
#' @md
#' @export
get_contour_plot <- function(img, min_val = 0.1){
  # min_val :
  img_contour_dt <- imager::imgradient(img, "xy") %>% imager::enorm %>%
    as.data.frame() %>% tidyverse::mutate(y_flip= -y)
  # imgradient : calculating the image gradients in 2D
  # enorm : distance
  # y_flip is only for reference purpose
  plot(img_contour_dt[img_contour_dt$value > min_val, 'x'],
       img_contour_dt[img_contour_dt$value > min_val, 'y_flip'],  ann = F)
}

#' edge filter x
#' @param dat
#' @md
#' @export
#' @md
#' @export

edge_filter_x <- function(dat){
  grid_dat <- dat[,'x'] %>% unique  # range of x
  minimum <- data.frame()
  maximum <- data.frame()
  for(grid_index in grid_dat){  # max/min points for each point in range of x
    minimum[grid_index, 'x'] <- grid_index
    maximum[grid_index, 'x'] <- grid_index
    minimum[grid_index, 'y'] <- min(dat[dat[,'x'] == grid_index, 'y'])
    maximum[grid_index, 'y'] <- max(dat[dat[,'x'] == grid_index, 'y'])
  }
  rbind(minimum, maximum)
}

#' edge filter y; same as edge_filter_y but for y
#' @param dat
#' @md
#' @export
edge_filter_y <- function(dat){
  grid_dat <- dat[,'y'] %>% unique
  minimum <- data.frame()
  maximum <- data.frame()
  for(grid_index in grid_dat){
    minimum[grid_index, 'x'] <- min(dat[dat[,'y'] == grid_index, 'x'])
    maximum[grid_index, 'x'] <- max(dat[dat[,'y'] == grid_index, 'x'])
    minimum[grid_index, 'y'] <- grid_index
    maximum[grid_index, 'y'] <- grid_index
  }
  rbind(minimum, maximum)
}

#' get the edge
#' edge filter y; same as edge_filter_y but for y
#' @param img image source
#' @param min_val threshold
#' @param type
#' @md
#' @export
get_edge <- function(img, min_val = 0.01, type){
  edge_list <- list()
  for(image_index in seq_along(img)){
    edge_dat <- imgradient(img[[image_index]], "xy") %>%
      enorm %>% as.data.frame() #  data.frame having the size of image variation
    edge_x <- edge_filter_x(edge_dat[edge_dat$value > min_val, c('x', 'y')]) # get the edge with image change larger than certain size (min_val) for x
    edge_y <- edge_filter_y(edge_dat[edge_dat$value > min_val, c('x', 'y')]) # for y
    edge_list[[paste0(type, '_edge_', image_index)]] <- na.omit(rbind(edge_x, edge_y)) # remove the NA's
  }
  edge_list
}

#' Eucledean distance
#' edge filter y; same as edge_filter_y but for y
#' @param x1
#' @param x2
#' @md
#' @export
euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

#' calculating the Eucledean distance from the center of the debris to the edge
#' @param x
#' @param y
#' @md
#' @export
#'
piece_dist <- function(x, y){
  distance <-c()
  for(i in 1:dim(y)[1]){
    distance <- c(distance, euc_dist(x, y[i,]))
  }
  distance
}

#' calcuating the distance from the center tot he edge for each sample
#' @md
#' @export
#'
get_piece_dist <- function(edge_list, type){
  dist_list <- list()
  for(image_index in seq_along(edge_list)){
    center <- c((max(edge_list[[image_index]]$x) + min(edge_list[[image_index]]$x)) / 2 ,
                (max(edge_list[[image_index]]$y) + min(edge_list[[image_index]]$y)) / 2) # center
    dist_list[[paste0(type, '_dist_', image_index)]] <- piece_dist(center, edge_list[[image_index]]) # distance from the center to the boundary
  }
  names(dist_list) <- str_to_lower(names(dist_list))
  dist_list
}

#' data frame for each sample
#' @md
#' @export
#'
get_dist_dt <- function(dist_list, type){
  dist_dt <- data.frame()
  for(image_index in seq_along(dist_list)){
    dist_dt <- rbind(dist_dt,
                     data.frame(
                       distance = dist_list[[image_index]],
                       key = type, sample_id = as.character(image_index)))
  }
  dist_dt
}



#' function to calculate the local curvature (EBImage vignettes)
#' @md
#' @export
#'
get_curvature <- function(img, h = 10){
  contours <-  ocontour(img)  # extract the image boundaries
  lc_out <- localCurvature(x=contours[[1]], h = h) # local curvature with smoothing window = h
  ## localCurvature() : contour, curvature, length
  lc_out$curvature # return curvature
}


#' function to assign the local curvature to the boundary for visulization
#' @md
#' @export
get_curvature_image <- function(img, h = 10){
  contours <-  ocontour(img)
  lc_out <- localCurvature(x=contours[[1]], h = h)
  index <-  lc_out $curvature >= 0
  pos <- neg <- array(0, dim(img))
  pos[lc_out$contour[index,]+1] <- lc_out$curvature[index]  # assign the curvature value to the boundary considering the signs of curvatures (positive/nevative)
  neg[lc_out$contour[!index,]+1] <- -lc_out$curvature[!index]
  out <- pos - neg
  out <- out %>% as_tibble %>%
    mutate(y_index = 1:nrow(out)) %>%
    gather(x_index, curvature, - y_index) %>%
    mutate(x_index = str_remove(x_index, pattern = "V")) %>%
    mutate(x_index = as.numeric(x_index)) %>%
    filter(abs(curvature) > 0) %>%
    mutate(curvature = ifelse(curvature >= 0, "Positive", "Negative") %>% as.factor)
  image_out <- ggplot(out) + aes(x = x_index, y = y_index, fill = curvature) + geom_tile() +
    theme_bw()
  image_out
}

#' restructure by samples
#' @md
#' @export
df_by_sample <- function(img_list, plastic_type){
  tbl_list <- list()
  for (image_index in seq_along(img_list)){
    tbl_list[[image_index]] <- tibble(curvature = img_list[[image_index]],
                                      sample_id = image_index,
                                      key = plastic_type,
                                      id = 1:length(img_list[[image_index]])
    )
  }
  bind_rows(tbl_list)
}

#' K is the length of the partition
#' @md
#' @export
perm_test <- function(sample_combined, K = 200){
  ref_tbl <-
    sample_combined %>%
    group_by(key, sample_id) %>%
    dplyr::summarise(max_id = max(id) - K + 1)

  sample_combined_max <- dplyr::left_join(sample_combined, ref_tbl, by = c("sample_id", "key"))
  candidates_pool <- sample_combined_max %>%
    filter(id <= max_id)
  picked_start_point <- candidates_pool %>%
    group_by(key, sample_id) %>%
    sample_n(1) %>%
    rename(start_point = id) %>%
    select(-curvature, -max_id)

  tbl_random_partition <-
    left_join(sample_combined, picked_start_point, by = c("sample_id", "key")) %>%
    mutate(end_point = start_point + K )
  out <- tbl_random_partition %>%
    filter(id >= start_point) %>%
    filter(id < end_point) %>%
    group_by(key, sample_id) %>%
    summarise(M = mean(abs(diff(curvature)))) %>%
    group_by(key) %>%
    summarise(M = mean(M))
  out
}
