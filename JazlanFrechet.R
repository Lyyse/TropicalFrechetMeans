library(combinat)


# Auxiliary functions


trop_dist <- function (x, y) {
  # Returns the tropical distance between points x and y.
  return (max(x - y) - min(x - y))
}


sum_dist <- function(ref, points, power=1) {
  # Returns the sum of tropical distances from ref to points.
  dists <- apply(points, 1, \(x) trop_dist(x, ref)**power)
  return(sum(dists))
}


sum_squares <- function(ref, points) {
  # Returns the sum of squared tropical distances from ref to points.
  squares <- apply(points, 1, \(x) trop_dist(x, ref)**2)
  return(sum(squares))
}


trop_ball <- function (centre, radius) {
  # Returns the vertices of a tropical ball.
  
  directions <- unique(permn(c(1,-1, rep(0, length(centre)-2))))
  
  return (lapply(directions, \(x) x*radius/2 + centre))
  
}


ssq_subgrad <- function (ref, points) {
  # Returns a subgradient at ref of the sum of squared tropical distances 
  # from ref to the given points.
  
  num_points <- dim(points)[1]
  dimension <- dim(points)[2]
  
  # Find all coordinate pairs that attain the maximum
  helper <- function (p) {
    i <- which.max(ref - p)
    j <- which.min(ref - p)
    return (2 * (ref[i] - ref[j] - p[i] + p[j]) *
              replace(numeric(dimension), c(i,j), c(1,-1)))
  }
  
  subgrads <- apply(points, 1, helper)
  subgrad <- rowSums(subgrads)
  return (subgrad)
  
}


# Fréchet mean computation


FM_hypersurface <- function (points) {
  # Returns the 1-skeleton of the FM hypersurface of a list of sample points.
  
  # Pass information to polymake
  points_str <- paste(apply(sample_points, 1, \(x) paste(x, collapse = ",")), 
                      collapse="] [")
  linux_command <- paste0("polymake --script test [", points_str, "]")
  result <- system2(command = "wsl", args = c(linux_command))
  
  # Process the vertices
  verts <- read.table("vertex_info.txt", sep=" ")
  colnames(verts)[1] <- "type"
  
  # Process the edges
  edges <- readLines("edge_info.txt")
  edges <- gsub("\\{|\\}", "", edges)  # Remove curly braces
  edges <- strsplit(edges, " ")      # Split by space
  edges <- lapply(edges, as.double)
  edges <- t(as.data.frame(edges)) + 1
  rownames(edges) <- NULL
  
  return (list(verts, edges))
}


fermat_weber_relaxed <- function (all_vertices, all_edges, target, points) {
  # Returns the relaxed Fermat-Weber set given the vertices and edges of
  # an FM hypersurface, a target (optimum + relaxation), and the sample points.
  
  # Check if there is a relaxed FW vertex on each edge
  edge_check <- function(edge) {
    
    if (all_vertices$type[edge[2]] == 0 & 
        all_vertices$sums[edge[1]] < target) {
      return (1)  # first index is a valid vertex; second is a ray
      
    } else if (all_vertices$type[edge[1]] == 0 & 
               all_vertices$sums[edge[2]] < target) {
      return (2)  # second index is a valid vertex; first is a ray
      
    } else if ((all_vertices$sums[edge[1]] < target & 
                all_vertices$sums[edge[2]] > target) |
               (all_vertices$sums[edge[2]] < target & 
                all_vertices$sums[edge[1]] > target)) {
      return (3)  # segment
    } else {
      return (0)  # there is no relaxed FW vertex here
    }
  }
  
  # Find a relaxed FW vertex on a valid segment
  segment_find_vert <- function(edge) {
    
    v1 <- all_vertices[edge[1],-c(1,2)]  # remove type and sums
    v2 <- all_vertices[edge[2],-c(1,2)]
    
    v1_sum <- all_vertices$sums[edge[1]]
    v2_sum <- all_vertices$sums[edge[2]]
    
    dir_vec <- v2 - v1
    proportion <- (target - v1_sum) / (v2_sum - v1_sum)
    valid_vert <- v1 + dir_vec * proportion
    
    return(valid_vert)
    
  }
  
  # Find a relaxed FW vertex on a valid ray
  ray_find_vert <- function(edge, class) {
    
    if (class == 1) {
      v0 <- all_vertices[edge[1],-c(1,2)]  # remove type and sums
      ray <- all_vertices[edge[2],-c(1,2)]
      v0_sum <- all_vertices$sum[edge[1]]
      
    } else if (class == 2) {
      v0 <- all_vertices[edge[2],-c(1,2)]
      ray <- all_vertices[edge[1],-c(1,2)]
      v0_sum <- all_vertices$sum[edge[2]]
      
    } else {
      stop("Improper class (only accepts 1 or 2)")
    }
    
    extra_sum <- sum_dist(v0 + ray, points)
    proportion <- (target - v0_sum) / (extra_sum - v0_sum)
    valid_vert <- v0 + ray * proportion
    
    return(valid_vert)
    
  }
  
  # Find a relaxed FW vertex on a valid ray
  check <- apply(all_edges, 1, \(x) edge_check(x))
  all_edges <- cbind(all_edges, check)
  
  FW_verts <- list()
  FW_verts <- rbind(FW_verts, all_vertices[all_vertices$sums == target])
  
  for (i in 1:nrow(all_edges)) {
    
    if (all_edges[i,3] == 1) {
      FW_verts <- rbind(FW_verts, ray_find_vert(all_edges[i,-3],1))
      
    } else if (all_edges[i,3] == 2) {
      FW_verts <- rbind(FW_verts, ray_find_vert(all_edges[i,-3],2))
      
    } else if (all_edges[i,3] == 3) {
      FW_verts <- rbind(FW_verts, segment_find_vert(all_edges[i,-3]))
      
    }
  }
  
  return (FW_verts)
  
}


fermat_weber_two_step <- function (points, delta=0.05) {
  # Returns the vertices of the relaxed FW polytope of a
  # list of sample points.
  
  # Initialise constants
  k <- nrow(points)
  
  # Step 0: no relaxation
  my_hyper <- FM_hypersurface(points)
  all_vertices <- my_hyper[[1]]  # i.e. pseudovertices
  all_edges <- my_hyper[[2]]     # Also includes rays
  
  sums <- apply(all_vertices, 1, 
                \(x) ifelse(x[1] == 1, sum_dist(x[-1], points), Inf))
  all_vertices <- cbind(sums, all_vertices)
  
  FW_optimum <- min(all_vertices$sums)
  vanilla_FW <- all_vertices[all_vertices$sums == FW_optimum, 
                             -ncol(all_vertices)]
  
  print("Step 0 done")
  
  # Step 1: consistent but sub-optimal relaxation
  relax1 <- FW_optimum * (log(k) / k) ** 0.5
  target1 <- FW_optimum + relax1
  relaxed_FW1 <- fermat_weber_relaxed(all_vertices, all_edges, 
                                      FW_optimum + relax1, points)
  
  # Compute the diameter of the relaxed FW set
  pairwise_dist <- matrix(NA, nrow = nrow(relaxed_FW1),
                              ncol = nrow(relaxed_FW1))
  for (i in 1:nrow(relaxed_FW1)) {
    for (j in i:nrow(relaxed_FW1)) {
      pairwise_dist[i,j] <- trop_dist(relaxed_FW1[i,], relaxed_FW1[j,])
    }
  }
  diameter <- max(pairwise_dist, na.rm = TRUE)
  
  print("Step 1 done")
  
  # Step 2: near-optimal relaxation
  relax2 <- (1 + delta) * diameter * (2 * log(log(k)) / k) ** 0.5
  target2 <- FW_optimum + relax2
  relaxed_FW2 <- fermat_weber_relaxed(all_vertices, all_edges, 
                                      FW_optimum + relax2, points)
  
  print("Step 2 done")
  
  return (relaxed_FW2)
  
}


greedy_frechet_mean <- function (points, perturb=0.001, power=2) {
  # Returns a tropical Fréchet mean of a list of sample points.
  # Uses Bo's descent method.
  # Specify perturb for the 'resolution'
  
  # Initialise constants
  num_points <- dim(points)[1]
  dimension <- dim(points)[2]
  perturb_vecs <- t(data.frame(unique(permn(c(1,-1, rep(0, dimension-2))))))
  
  # Begin at the specified starting point
  frechet_mean <- points[1,]
  sod <- sum_dist(frechet_mean, points, power)
  
  while (TRUE) {
    # Perturb the current point and compute the sum of squares
    # Inefficient, recomputes a lot of the same points
    perturbed <- apply(perturb_vecs, 1, \(x) x + frechet_mean)
    new_sod <- apply(perturbed, 2, \(x) sum_dist(x, points, power))
    
    # Check if the sum of squares can be reduced
    if (all(sod <= new_sod)) {
      # Saddle point -> how to decide which direction to go?
      break
    }
    
    # Move as far as we can in the direction of the decrease.
    direction <- as.numeric(perturb_vecs[which.min(new_sod),])
    while (TRUE) {
      new <- frechet_mean + direction
      new_sod <- sum_dist(new, points, power)
      if (sod <= new_sod) {
        break
      }
      frechet_mean <- new
      sod <- new_sod
    }
  }
  
  return (list(frechet_mean, sod))
}


subgrad_frechet_mean <- function (points) {
  # Returns a tropical Fréchet mean of a list of sample points.
  # Uses the subgradient method.
  
  # Initialise constants
  num_points <- dim(points)[1]
  dimension <- dim(points)[2]
  frechet_mean <- rep(0, dimension)
  count <- 1
  
  # Find the subgradient of the current point
  while (TRUE) {
    current_sos <- sum_squares(frechet_mean, points)
    subgrad <- ssq_subgrad(frechet_mean, points)
    new <- frechet_mean - subgrad / count
    print(sum_squares(new, points))
    if (new == frechet_mean) {
      break
    } else {
      frechet_mean = new
    }
    count <- count + 1
    
    # TODO stopping condition, more thoughtful choice of subgradient?
  }
}


frechet_polytope <- function (points) {
  # Returns a matrix whose rows are the vertices of the Fréchet mean polytope 
  # of a list of sample points.
  
  # Initialise constants
  num_points <- dim(points)[1]
  dimension <- dim(points)[2]
  
  # Find one Fréchet mean of the sample
  frechet_mean <- greedy_frechet_mean(points)[[1]]
  
  # Construct tropical balls around each point
  radius <- trop_dist(frechet_mean, points[1,])
  
  A_mat <- 0  # This leaves some redundant inequalities
  b_vec <- 0
  
  for (k in 1:num_points) {
    ball <- trop_ball(points[k,], radius)
    A_mat <- rbind(A_mat, ball[[1]])
    b_vec <- c(b_vec, ball[[2]])
  }
  
  # Find the vertices of the intersection of these balls
  # vertexenum is not available in my version of R!!!!
  # TODO
}


# Tests

sample1 <- c(0,6,1,7,9,10)
sample2 <- c(0,0,0,8,2,3)
sample3 <- c(0,3,-5,9,10,11)
sample4 <- c(0,3,6,1,5,6)
sample_points <- rbind(sample1, sample2, sample3, sample4)
greedy_frechet_mean(sample_points, perturb=0.001, 1)

relaxed_plot_verts <- fermat_weber_two_step(sample_points)
relaxed_plot_verts
