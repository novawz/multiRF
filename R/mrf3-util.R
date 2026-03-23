## Utility functions for mrf3.R

# Pre-process nativeArray once and split by treeID
.prep_tree_dfs <- function(mod) {
  # Native C++ engine: build nativeArray-equivalent from tree_info
  if (!is.null(mod$tree_info) && is.null(mod$forest$nativeArray)) {
    return(.prep_tree_dfs_native(mod))
  }

  na_raw <- mod$forest$nativeArray
  na_df <- data.frame(na_raw)
  node.stat <- mod$node.stats
  if (!is.null(node.stat) && NROW(node.stat) == nrow(na_df)) {
    for (col in setdiff(colnames(node.stat), colnames(na_df))) {
      na_df[[col]] <- node.stat[[col]]
    }
  }
  if (!"nodeSZ" %in% names(na_df)) na_df$nodeSZ <- 1L
  if (!"dpthST" %in% names(na_df)) {
    # Compute depth proxy when node.stats unavailable (rfsrc >= 3.5)
    na_df$dpthST <- ave(na_df$nodeID, na_df$treeID,
                        FUN = function(ids) seq_along(ids))
  }
  # Split once — returns a named list of per-tree data frames
  split(na_df, na_df$treeID)
}

# Convert native engine tree_info to nativeArray-like per-tree data frames.
# The C++ engine stores nodes in BFS order, but build_tree_network_cpp
# expects pre-order (DFS) traversal, so we reorder here.
.prep_tree_dfs_native <- function(mod) {
  tree_info <- mod$tree_info
  nt <- length(tree_info)
  dfs <- vector("list", nt)
  for (t in seq_len(nt)) {
    ti <- tree_info[[t]]
    n_nodes <- length(ti$split_var)
    # Compute pre-order traversal of the tree
    dfs_order <- .bfs_to_preorder(ti$left, ti$right, n_nodes)
    # split_var: 0-indexed X column (-1 = leaf) -> parmID: 1-indexed (0 = leaf)
    parmID <- ifelse(ti$split_var[dfs_order] < 0, 0L, ti$split_var[dfs_order] + 1L)
    dfs[[t]] <- data.frame(
      treeID = rep(t, n_nodes),
      nodeID = seq_len(n_nodes),
      parmID = parmID,
      nodeSZ = as.integer(ti$nodesize[dfs_order]),
      dpthST = as.integer(ti$depth[dfs_order]),
      stringsAsFactors = FALSE
    )
  }
  names(dfs) <- as.character(seq_len(nt))
  dfs
}

# Convert BFS-ordered tree to pre-order (DFS) traversal indices.
# left/right are 0-indexed child indices (-1 = no child / leaf).
.bfs_to_preorder <- function(left, right, n_nodes) {
  order_vec <- integer(n_nodes)
  stack <- 1L  # 1-indexed root
  idx <- 0L
  while (length(stack) > 0L) {
    node <- stack[length(stack)]
    stack <- stack[-length(stack)]
    idx <- idx + 1L
    order_vec[idx] <- node
    # Push right first so left is processed first (stack = LIFO)
    r <- right[node]
    l <- left[node]
    if (r >= 0) stack <- c(stack, r + 1L)  # 0-indexed -> 1-indexed
    if (l >= 0) stack <- c(stack, l + 1L)
  }
  order_vec
}

# Get tree net from random forest
get_tree_net <- function(mod, tree.id, tree_dfs = NULL){

  xvar.names <- mod$xvar.names
  xvar.factor <- mod$xvar.factor

  if (!is.null(tree_dfs)) {
    tree.df <- tree_dfs[[as.character(tree.id)]]
  } else {
    # Fallback: compute on the fly (slow, for backwards compat)
    if (!is.null(mod$tree_info) && is.null(mod$forest$nativeArray)) {
      # Native engine: build single-tree df from tree_info (DFS order)
      ti <- mod$tree_info[[tree.id]]
      n_nodes <- length(ti$split_var)
      dfs_order <- .bfs_to_preorder(ti$left, ti$right, n_nodes)
      parmID <- ifelse(ti$split_var[dfs_order] < 0, 0L, ti$split_var[dfs_order] + 1L)
      tree.df <- data.frame(
        treeID = rep(tree.id, n_nodes),
        nodeID = seq_len(n_nodes),
        parmID = parmID,
        nodeSZ = as.integer(ti$nodesize[dfs_order]),
        dpthST = as.integer(ti$depth[dfs_order]),
        stringsAsFactors = FALSE
      )
    } else {
      na_df <- data.frame(mod$forest$nativeArray)
      node.stat <- mod$node.stats
      if (!is.null(node.stat) && NROW(node.stat) == nrow(na_df)) {
        for (col in setdiff(colnames(node.stat), colnames(na_df))) {
          na_df[[col]] <- node.stat[[col]]
        }
      }
      if (!"nodeSZ" %in% names(na_df)) na_df$nodeSZ <- 1L
      if (!"dpthST" %in% names(na_df)) na_df$dpthST <- 1L
      tree.df <- na_df[na_df$treeID == tree.id, , drop = FALSE]
    }
  }

  #converted.tree <- display.tree
  vars.id <- data.frame(var = c("<leaf>", xvar.names), parmID = 0:length(xvar.names), stringsAsFactors = FALSE)
  tree.df$var <- vars.id$var[match(tree.df$parmID, vars.id$parmID)]

  num.node <- tree.df$nodeID %>% unique %>% length

  var.count <- integer(nrow(tree.df))
  for (v in unique(tree.df$var)) {
    pt <- tree.df$var == v
    var.count[which(pt)] <- seq_len(sum(pt))
  }

  tree.df$var_count <- var.count
  tree.df$var_conc <- paste0(tree.df$var, "_", tree.df$var_count)

  n_internal <- sum(tree.df$var != "<leaf>")
  n_leaf     <- sum(tree.df$var == "<leaf>")
  tree.df$var.tip.id <- integer(nrow(tree.df))
  tree.df$var.tip.id[tree.df$var != "<leaf>"] <- n_leaf + seq_len(n_internal)
  tree.df$var.tip.id[tree.df$var == "<leaf>"] <- tree.df$var_count[tree.df$var == "<leaf>"]

  is_leaf_vec <- tree.df$var == "<leaf>"

  network <- build_tree_network_cpp(
    var_conc   = tree.df$var_conc,
    var_tip_id = tree.df$var.tip.id,
    nodeSZ     = tree.df$nodeSZ,
    dpthST     = tree.df$dpthST,
    is_leaf    = is_leaf_vec
  )

  network <- data.frame(network, stringsAsFactors = FALSE)

  return(network)

}

# Get leaf node information
get_leaf_ds <- function(mod, tree.membership, net){
  
  # Find upper level node
  prev_leaf <- dplyr::filter(.data = net, is_leaf == 1)
  
  # Update mem_id and membership
  net$mem_id_old <- net$mem_id
  mem_update <- slice_max(.data = prev_leaf, order_by = mem_id, by = from, with_ties = FALSE)
  mem_old <- slice_min(.data = prev_leaf, order_by = mem_id, by = from, with_ties = FALSE)
  mem_new <- setNames(mem_update$mem_id, mem_update$from)
  mem_org_id <- mem_old$mem_id
  names(mem_org_id) <- mem_new
  sapply(1:nrow(mem_old), function(id) net$mem_id[net$mem_id == mem_org_id[id]] <<- as.numeric(names(mem_org_id)[id]))
  
  find_node <- unique(prev_leaf$from)
  
  leaf_select <- filter(net, to %in% find_node)
  map_id <- unique(leaf_select$from)
  
  leaf_select <- filter(net, from %in% map_id)
  leaf_select$mem_id[match(names(mem_new), leaf_select$to)] <- mem_new
  drop_leaf <- leaf_select$from[leaf_select$mem_id == 0]
  leaf_select <- leaf_select[!leaf_select$from %in% drop_leaf,]
  
  # Update leaf information
  
  net$is_leaf[net$to %in% prev_leaf$to] <- 0
  net$is_leaf[net$from %in% leaf_select$from] <- 1
  
  net$mem_id[match(names(mem_new), net$to)] <- mem_new
  
  tree_mem <- tree.membership
  
  sapply(1:length(mem_org_id),
         function(id){
           tree_mem[tree_mem == mem_org_id[id]] <<-
             as.numeric(names(mem_org_id)[id])
         })
  
  
  if(any(leaf_select$mem_id %in% 0)) {
    from_id <- leaf_select[leaf_select$mem_id %in% 0,"from"]
    leaf_select[leaf_select$from %in% from_id,"is_leaf"] <- 0
    leaf_drop <- leaf_select[leaf_select$is_leaf == 0,]
    net[net$to %in% leaf_drop$to,"is_leaf"] <- 0
    net[net$to %in% leaf_drop$to,"leaf_stack"] <- 1
  }
  
  return(
    list(net = net,
         tree.mem = tree_mem)
  )
}

# Get response outcome splitting scores for splits in a tree
get_Y_imp <- function(net, tree.membership, dat,
                      robust = FALSE, w = NULL, ytry = 1, seed = -5,
                      top = 5){

  node_id <- unique(net$from)
  dat_mat <- as.matrix(dat)
  mem_int <- as.integer(tree.membership)

  var_imp_ls <- lapply(node_id, function(id) {
    children <- net[net$from == id, , drop = FALSE]
    mem_id <- children$mem_id

    # Use C++ for the heavy split statistic computation
    split_stat <- compute_split_stats_cpp(
      Y = dat_mat,
      membership = mem_int,
      mem_left  = as.integer(mem_id[1]),
      mem_right = as.integer(mem_id[2])
    )

    if (robust) {
      split_stat_raw <- split_stat
      keep_idx <- order(split_stat, decreasing = TRUE)[seq_len(min(top, length(split_stat)))]
      mask <- rep(0, length(split_stat))
      mask[keep_idx] <- 1
      split_stat <- split_stat * mask
      if (!any(is.finite(split_stat) & split_stat != 0)) {
        split_stat <- split_stat_raw
      }
    }

    idx <- which.max(split_stat)
    varY <- mean(split_stat)
    varY_all <- split_stat
    names(varY_all) <- colnames(dat)
    names(varY) <- colnames(dat)[idx]

    list(varY = varY, varY_all = varY_all)
  })

  var_imp <- unlist(var_imp_ls %>% purrr::map("varY"))
  var_imp_all <- Reduce(rbind, var_imp_ls %>% purrr::map("varY_all"))
  
  list(var_imp = var_imp, var_imp_all = var_imp_all)
  
}

# Get importance for leaves in a tree
get_tree_imp <- function(mod, dat = NULL, robust = FALSE, tree.membership, net, calc = "Both", M = NULL, w = NULL, ytry = 1, weighted = FALSE, seed = -5){

  if(is.null(dat)){
    dat <- mod$xvar
  }

  if(calc %in% c("Y", "Both")){
    datY <- mod$yvar[rownames(dat),]
  }

  updated_net <- net

  top_node_info <- filter(.data = updated_net, is_leaf == 1)

  node_ds <- setNames(top_node_info$inv_d, top_node_info$from)

  old_net <- updated_net[updated_net$from %in% top_node_info$from,]
  scores_imp <- NULL

  old_net_corr <- group_by(.data = old_net, from)
  old_net_corr <- dplyr::summarise_at(old_net_corr, .vars = "inv_d", .funs = mean)

  node_ds <- node_ds[old_net_corr$from]

  use_sample <- old_net_corr$from

  match_old_net <- old_net[old_net$from %in% use_sample,]

  # Get important variable
  # Get Y

  if(calc %in% c("Both", "Y")){

    impY_ls <- get_Y_imp(net = match_old_net, tree.membership = tree.membership, robust = robust, dat = datY, ytry = ytry, seed = seed)
    impY <- impY_ls$var_imp
    updated_net$Y_id[unique(match(match_old_net$from, updated_net$from))] <- names(impY)
    scores_impY <- updated_net$inv_d[match(unique(match_old_net$from), updated_net$from)]
    if(weighted) {
      scores_impY <- scores_impY * updated_net$edge[match(unique(match_old_net$from), updated_net$from)]
    }
    names(scores_impY) <- names(impY)

    impY_mat <- impY_ls$var_imp_all * impY
    # scores_impY <- scores_impY
  }

  # Get X
  if(calc %in% c("Both", "X")){

    imp_var <- top_node_info[match(use_sample, top_node_info$from),"from"]
    if(robust) {
      impX_ls <- get_Y_imp(net = match_old_net, tree.membership = tree.membership, robust = robust, dat = dat, ytry = ytry, seed = seed)
      impX <- impX_ls$var_imp
    }

    sub_net <- old_net_corr[match(unique(match_old_net$from), old_net_corr$from),]
    scores_imp <- (sub_net$inv_d)
    if(weighted) {
      scores_imp <- scores_imp * old_net_corr$edge
    }

    # scores_imp <- drop_case[match(imp_var,drop_case$to),] %>% pull(corr)
    names(scores_imp) <- gsub("^(.*)_.*", "\\1",imp_var[match(sub_net$from, imp_var)])
    if (robust) impX_mat <- impX_ls$var_imp_all * (sub_net$inv_d) * impX
  }
  if (robust) {
    if(!is.null(impY) | !is.null(impX)) {
      
      if(is.null(dim(impY_mat))) {
        M <- M + sub_net$inv_d * tcrossprod(impX_mat, impY_mat)
      } else {
        M_s <- lapply(1:nrow(impY_mat), function(k) {
          x <- impX_mat[k,]
          y <- impY_mat[k,]
          sub_net$inv_d[k] * tcrossprod(x,y)
        })
        M <- M + Reduce("+", M_s)
      }
    }
    #M[names(scores_imp),] <- M[names(scores_imp),] + impY_mat
    #M[,names(scores_impY)] <- M[,names(scores_impY)] + t(impX_mat)
  }  else { M <- NULL }
  
  if(calc == "Both"){
    scores_imp <- list(X = scores_imp, Y = scores_impY)
  } else if(calc == "X"){
    scores_imp <- list(X = scores_imp)
  } else {
    scores_imp <- list(Y = scores_impY)
  }


  return(
    list(
      net = updated_net,
      dat = dat,
      imp_var = scores_imp,
      M = M
    )
  )
}

# Update tree importance from bottom to top
update_iter_imp <- function(mod, tree.id, calc = "Both", robust = FALSE, w = NULL, ytry = 1, weighted = FALSE, seed = -5,
                            tree_dfs  = NULL) {

  # Get the tree structure for the specified tree.id
  net <- get_tree_net(mod = mod, tree.id = tree.id, tree_dfs = tree_dfs)

  # Get the membership information for the specified tree.id
  mem <- mod$membership[, tree.id]

  # Initialize an empty data frame 'dat'
  dat <- NULL

  # Initialize an index 'i' and k
  i <- 1
  k <- 1

  # Initialize an empty list 'imp'
  imp <- list()

  net$Y_id <- NA
  if(robust) {
    mat <- matrix(0, ncol = ncol(mod$yvar), nrow = ncol(mod$xvar))
    
    colnames(mat) <- colnames(mod$yvar)
    rownames(mat) <- colnames(mod$xvar)
  }
  else mat <- NULL 


  # Iterate until all nodes in the tree are leaves or there are only 2 or fewer nodes left
  while (any(is.null(net$is_leaf), sum(net$used_leaf == 0) > 0, is.null(net$used_leaf))) {

    if(k > 1){
      if(is.null(net$used_leaf)) {
        net$used_leaf <- ifelse(net$is_leaf == 1, 1, 0)
      } else{
        net$used_leaf[net$is_leaf == 1] <- 1
      }
      # Get initial leaf information
      update_net <- get_leaf_ds(mod = mod,
                        tree.membership = mem,
                        net = net)
      mem <- update_net$tree.mem
      net <- update_net$net
    } else {
      # Create new columns
      if(is.null(net$is_leaf)) {
        net$is_leaf <- ifelse(grepl("leaf" ,net$to), 1, 0)
      }
      if(is.null(net$mem_id)) {
        net$mem_id <- ifelse(net$is_leaf == 1, net$to_id, 0)
      }
      find_node <- net[net$is_leaf == 1,]
      find_node <- dplyr::summarise(.data = find_node, n = n(), .by = from)
      find_node <- filter(.data = find_node, n == 1)
      net[net$is_leaf == 1,"is_leaf"][net[net$is_leaf == 1,]$from %in% find_node$from] <- 0
    }

    # Get the updated tree importance statistics
    update_ls <- get_tree_imp(
      mod = mod,
      dat = dat,
      robust = robust,
      tree.membership = mem,
      net = net,
      calc = calc,
      M = mat,
      w = w,
      ytry = ytry,
      weighted = weighted,
      seed = seed
    )

    # Update 'net', 'mem', and 'dat' with the results from 'update_ls'
    net <- update_ls$net
    dat <- update_ls$dat
    
    mat <- update_ls$M

    # If importance for a variable is available, add it to the 'imp' list
    if (!is.null(update_ls$imp_var)) {
      imp[[i]] <- update_ls$imp_var
      # if(i == 1) {
      #   varls <- update_ls$imp_var
      #   varls <- purrr::map(varls, ~.-.)
      #   imp[[i]] <- varls
      # }
      i <- i + 1
    }

    k <- k + 1

  }


  # Filter out the root node(s)
  root_net <- filter(net, is_leaf == 1)

  # Initialize variables for variable importance for X and Y
  # if (calc %in% c("X", "Both")) {
  #
  #   imp_root <- 1/2
  #   names(imp_root) <- gsub("^(.*)_.*", "\\1", unique(root_net$from))
    impX <- unlist(purrr::map(imp, "X"))
  #   impX <- c(impX, imp_root)
  # }
  #
  if (calc %in% c("Y", "Both")) {
    impY <- unlist( purrr::map(imp, "Y") )
  }

  # Create lists of variable importance and variable names
  if (calc == "Both") {
    imp_ls <- list(impX, impY)
    var_name <- list(colnames(mod$xvar), colnames(mod$yvar))
  } else if (calc == "X") {
    imp_ls <- list(impX)
    var_name <- list(colnames(mod$xvar))
  } else {
    imp_ls <- list(impY)
    var_name <- list(colnames(mod$yvar))
  }

  # Calculate variable importance scores for each variable
  imp_col <- plyr::llply(
    1:length(imp_ls),
    .fun = function(l) {
      get_iv(var_name[[l]], imp_ls[[l]])
    }
  )

  # Set appropriate names for the variable importance scores
  if (calc == "Both") {
    names(imp_col) <- c("X", "Y")
  } else if (calc == "X") {
    names(imp_col) <- "X"
  } else {
    names(imp_col) <- "Y"
  }

  # add penalty
  # x_freq <- mod$var.used[tree.id,]
  # x_freq <- x_freq[x_freq != 0]
  # imp_col <- add_lambda(imp_col, net, x_freq, lambda)

  out <- list(imp_ls = imp_col, net = net, M = mat)
  return(out)
}

# Add penalty to weights
add_lambda <- function(imp_ls, net, x_freq, lambda){

  if(all(is.na(net$Y_id))){
    freq <- list(X = x_freq)
  } else {
    freq <- list(X = x_freq, Y = table(net$Y_id))
  }

  new_ww <- plyr::llply(
    names(imp_ls),
    .fun = function(i){

      var <- freq[[i]]
      var_imp_all <- imp_ls[[i]]
      var_imp <- var_imp_all[var_imp_all != 0]
      var <- var[names(var_imp)]

      if(any(var == 1)){
        var_imp[var == 1] <- var_imp[var == 1] * lambda
        var_imp_all[names(var_imp)] <- var_imp
      }

      var_imp_all
    }
  )

  names(new_ww) <- names(imp_ls)

  new_ww

}

#' Get forest importance
#' @param mod A fitted forest model from the native multiRF engine or the
#' optional `randomForestSRC` fallback.
#' @param parallel Logical; whether to parallelize across trees.
#' @param robust Logical; whether to use robust matrix-based aggregation.
#' @param calc Which importance side to compute: `"X"`, `"Y"`, or `"Both"`.
#' @param weighted Logical; whether to use weighted importance updates.
#' @param use_depth Logical; whether to aggregate non-zero depths instead of simple mean.
#' @param normalized Logical; whether to l2-normalize returned importance.
#' @param w Optional case weights.
#' @param ytry Response sampling proportion used in node-level updates.
#' @param cores Number of CPU cores used when `parallel = TRUE`.
#' @param seed Random seed passed to stochastic components.
#' @rdname get_imp_forest
get_imp_forest <- function(mod, parallel = FALSE, robust = FALSE, calc = "Both", weighted = FALSE, use_depth = FALSE, normalized = FALSE,
                           w = NULL, ytry = 1, cores = NULL, seed = -5){

  nt <- mod$ntree

  if(parallel){
    cores <- sanitize_mc_cores(cores = cores, fallback = 1L)
    if(Sys.info()["sysname"] == "Windows"){
      cluster <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cluster)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
    } else {
      doParallel::registerDoParallel(cores)
    }
  }

    if(calc == "Both") {
      cc <- "Both"
      idx <- c("X", "Y")
    }
    if(calc == "X") {cc <- "X";idx <- c("X")}
    if(calc == "Y") {cc <- "Y";idx <- c("Y")}

    # Pre-process nativeArray ONCE, split by tree — avoids per-tree overhead
    tree_dfs <- .prep_tree_dfs(mod)

    results <- plyr::llply(1:nt,
                       .fun = function(t){
                         update_iter_imp(mod,
                                         tree.id = t,
                                         calc = cc,
                                         robust = robust,
                                         ytry = ytry,
                                         w = w,
                                         weighted = weighted,
                                         seed = seed,
                                         tree_dfs  = tree_dfs)
                       }, .parallel = parallel)

    imp <- purrr::map(results, "imp_ls")
    net <- purrr::map(results, "net")

    if (robust) {

      M <- purrr::map(results, "M")
      imp <- Reduce("+", M)/nt
      impY <- apply(imp, 2, function(x) sqrt(sum(x^2))/nrow(imp)); impX <- apply(imp, 1, function(x) sqrt(sum(x^2))/ncol(imp))
      # M1 <- M1/max(M1); M2 <- M2/max(M2)
      # impX <- diag(M1); impY <- diag(M2)
      # impX <- rowMeans(M); impY <- colMeans(M)
      if(normalized){
        imp_ls <- list(X = impX/sqrt(sum(impX^2)), Y = impY/sum(sqrt(impY^2)))
      } else {
        imp_ls <- list(X = impX, Y = impY)
      }
      
    } else {
      imp_ls <- plyr::llply(idx,
                            .fun = function(im){
                              raw <- purrr::map(imp, im)
                              lens <- vapply(raw, length, integer(1))
                              raw <- raw[lens > 0L]
                              if (length(raw) == 0L) return(NULL)
                              im_df <- Reduce(cbind, raw)
                              if(!is.null(ncol(im_df))){
                                if(use_depth) {
                                  rs <- rowSums(im_df != 0)
                                  rs[rs == 0] <- -999
                                  imp <- rowSums(im_df)/rs
                                } else {
                                  imp <- rowMeans(im_df)
                                }
                              } else imp <- im_df
                              if(normalized){
                                imp <- imp/sqrt(sum(imp^2))
                              }
                              return(imp)
                            })
    }
    
    
  
    names(imp_ls) <- idx


  out <- list(
    imp_ls = imp_ls,
    imp_ls_init = imp,
    net = net
  )
  return(out)
}

# Match importance name to column name of the data frame
get_iv <- function(var_name, imp){

  var_v <- rep(0, length(var_name))
  names(var_v) <- var_name

  if (length(imp) == 0L || is.null(names(imp))) {
    return(var_v)
  }

  nm <- names(imp)
  ok <- !is.na(nm) & nzchar(nm) & is.finite(imp)
  if (!any(ok)) {
    return(var_v)
  }

  imp <- as.numeric(imp[ok])
  nm <- nm[ok]
  keep <- nm %in% var_name
  if (!any(keep)) {
    return(var_v)
  }

  imp_max <- tapply(imp[keep], nm[keep], max)
  var_v[names(imp_max)] <- as.numeric(imp_max)

  return(var_v)

}

#' Get multi-omics weights
#' @param mod_list A named list of fitted random forest models.
#' @param dat.list A named list of omics datasets used to align weights.
#' @param y Optional response vector (reserved for extensions).
#' @param weighted Logical; whether to use weighted importance updates.
#' @param use_depth Logical; whether to average depth-aware importances.
#' @param robust Logical; whether to use robust matrix-based aggregation.
#' @param parallel Logical; whether to parallelize across models.
#' @param normalized Logical; whether to normalize the merged weights.
#' @param calc Which importance side to compute: `"X"`, `"Y"`, or `"Both"`.
#' @param ytry Response sampling proportion used in node-level updates.
#' @param w Optional case weights.
#' @param cores Number of CPU cores used when `parallel = TRUE`.
#' @param seed Random seed passed to stochastic components.
#' @param ... Additional arguments for downstream helper functions.
#' @rdname get_multi_weights
get_multi_weights <- function(mod_list, dat.list, y = NULL, weighted = FALSE,  use_depth = FALSE, robust = FALSE,
                              parallel = FALSE, normalized = TRUE, calc = "Both", ytry = 1,
                              w = NULL, cores = NULL, seed = -5, ...){

  mod_names <- names(mod_list)

  # Fast path: if ALL models have pre-computed imd_weights (native engine),
  # skip the expensive post-hoc tree traversal entirely.
  all_have_imd <- all(vapply(mod_list, function(m) !is.null(m$imd_weights), logical(1)))
  if (all_have_imd) {
    message("  Using pre-computed IMD weights from native engine (zero extra cost).")
    weight_l <- lapply(mod_names, function(m_name) {
      iw <- mod_list[[m_name]]$imd_weights
      m_name_sep <- rev(unlist(stringr::str_split(m_name, "_")))
      names(iw) <- m_name_sep
      iw
    })

    block_names <- names(dat.list)
    weight_list <- lapply(block_names, function(bn) {
      ww <- purrr::compact(purrr::map(weight_l, bn))
      if (length(ww) == 0L) return(setNames(numeric(ncol(dat.list[[bn]])), colnames(dat.list[[bn]])))
      w_out <- Reduce("+", ww) / length(ww)
      if (normalized) {
        denom <- sqrt(sum(w_out^2))
        if (is.finite(denom) && denom > 0) w_out <- w_out / denom
      }
      w_out
    })
    names(weight_list) <- block_names

    # Build per-tree weight distributions (for method="test" in mrf3_vs)
    # Format: list of connections, each = list of ntree elements,
    #   each = named list(block1 = vec, block2 = vec)
    all_have_pt <- all(vapply(mod_list, function(m) !is.null(m$imd_weights_per_tree), logical(1)))
    weight_list_init <- NULL
    if (all_have_pt) {
      weight_list_init <- lapply(mod_names, function(m_name) {
        ipt <- mod_list[[m_name]]$imd_weights_per_tree  # list(X=mat, Y=mat)
        m_name_sep <- rev(unlist(stringr::str_split(m_name, "_")))
        nt_local <- ncol(ipt$X)
        # list of ntree elements, each = list(block1 = vec, block2 = vec)
        lapply(seq_len(nt_local), function(t) {
          out <- list(ipt$X[, t], ipt$Y[, t])
          names(out) <- m_name_sep
          out
        })
      })
    }

    return(list(
      weight_list = weight_list,
      weight_list_init = weight_list_init,
      net = NULL
    ))
  }

  # Slow path: post-hoc tree traversal (used when models come from rfsrc or
  # native engine without pre-computed IMD)
  results <- get_results(mod_list = mod_list, parallel = parallel, robust = robust, weighted = weighted, normalized = FALSE, use_depth = use_depth,
                         calc = calc, w = w, cores = cores, ytry = ytry, seed = seed)

  net <- purrr::map(results, "net")
  weight_l <- purrr::map(results, "wl")
  weight_l_init <- purrr::map(results, "wl_init")
  # M_list <- purrr::map(results, "M")
  
  if(length(weight_l) > 1){
    weight_list <- plyr::llply(
      names(dat.list),
      .fun = function(i){
        w <- purrr::map(weight_l, i)
        w <- purrr::compact(w)
        w <- (Reduce("+", w))/length(w)
        if(normalized) {
          denom <- sqrt(sum(w^2))
          if (is.finite(denom) && denom > 0) {
            w <- w/denom
          }
        }
        w
      }
    )
    names(weight_list) <- names(dat.list)
  } else {
    weight_list <- weight_l[[1]]
    wl <- plyr::llply(
      names(weight_list),
      .fun = function(i){
        w <- weight_list[[i]]
        if(normalized) {
          denom <- sqrt(sum(w^2))
          if (is.finite(denom) && denom > 0) {
            w <- w/denom
          }
        }
        w
      }
    )
    names(wl) <- names(weight_list)
    weight_list <- wl
    weight_list <- weight_list[names(dat.list)]
  }

  out <- list(weight_list = weight_list,
              weight_list_init = weight_l_init,
              net = net)

  return(out)
}


cal_freq <- function(mod, net){

  x_freq <- colMeans(mod$var.used != 0)
  freq <- list(X = x_freq)

  if(is.null(mod$yvar) || identical(class(mod)[3], "class+")){
    return(freq)
  } else {

    f <- rep(0, ncol(mod$yvar))
    names(f) <- colnames(mod$yvar)
    y_freq <- purrr::map(
      net,
      ~{
        id <- .[["Y_id"]]
        tb <- table(id)
        names(tb)
      }
    )
    y_freq <- table(unlist(y_freq))/mod$ntree
    f[names(y_freq)] <- y_freq

    freq <- c(freq,list(Y = f))

    return(freq)
  }


}

# step_two_weight <- function(two_step, rm_noise, normalized, weight_l, dat.list, freq_ls, s){
# 
#   if(two_step){
#     if(length(freq_ls) > 1){
#       freq <- plyr::llply(
#         names(dat.list),
#         .fun = function(i){
# 
#           fr <- purrr::map(freq_ls, i)
#           fr <- purrr::compact(fr)
#           Reduce("+", fr)/length(fr)
#         }
#       )
#       names(freq) <- names(dat.list)
#     } else {
#       freq <- freq_ls[[1]]
#     }
#   }
# 
#   if(length(weight_l) > 1){
#     weight_list <- plyr::llply(
#       names(dat.list),
#       .fun = function(i){
#         w <- purrr::map(weight_l, i)
#         w <- purrr::compact(w)
#         w <- (Reduce("+", w))/length(w)
#         if(two_step){
#           x <- freq[[i]]
#           # o <- sort(x, decreasing = T)[1:ceiling(a[i] * length(x))]
#           # x[x < min(o)] <- 0
#           w <- w * x
#         }
#         if(rm_noise){
#           x <- s * sd(w)
#           w[w < x] <- 0
#           w[w > x] <- w[w > x] - x
#           
#         } 
#         if(normalized) {
#           w <- w/sqrt(sum(w^2))
#         }
#         
#         w
#       }
#     )
#     names(weight_list) <- names(dat.list)
#   } else {
#     weight_list <- weight_l[[1]]
#     wl <- plyr::llply(
#       names(weight_list),
#       .fun = function(i){
#         w <- weight_list[[i]]
#         if(two_step){
# 
#           x <- freq[[i]]
#           # o <- sort(x, decreasing = T)[1:ceiling(a[i] * length(x))]
#           # x[x < min(o)] <- 0
#           w <- w * x
#         }
#         if(rm_noise){
#           x <- s * sd(w)
#           w[w < x] <- 0
#           w[w > x] <- w[w > x] - x
# 
#         } 
#         if(normalized) {
#           w <- w/sqrt(sum(w^2))
#         }
#         
#         w
#       }
#     )
#     names(wl) <- names(weight_list)
#     weight_list <- wl
#     weight_list <- weight_list[names(dat.list)]
#   }
# 
#   weight_list
# }

get_results <- function(mod_list, parallel,
                        normalized = FALSE, weighted = FALSE, robust = FALSE,
                        use_depth = FALSE, calc, w = NULL, ytry = 1, cores = NULL, seed = -5){

  mod_names <- names(mod_list)
  plyr::llply(
    mod_names,
    .fun = function(m_name){

      mod <- mod_list[[m_name]]
      # l <- lambda[m_name]
      if(!is.null(w)) {
        w0 <- w[[gsub("_.*", "", m_name)]]
      } else {w0 <- NULL}
      results <- get_imp_forest(mod, parallel = parallel, robust = robust, normalized = normalized,
                                weighted = weighted, calc = calc,  w = w0,
                                ytry = ytry, cores = cores, use_depth = use_depth, seed = seed)

      wl <- results$imp_ls
      wl_init <- results$imp_ls_init
      
      
      net <- results$net

      m_name_sep <- unlist(stringr::str_split(m_name, "_"))

      names(wl) <- rev(m_name_sep)
      # freq <- cal_freq(mod, net)
      # names(freq) <- rev(m_name_sep)
      return(list(
        wl = wl,
        wl_init = wl_init,
        net = net
        # freq = freq
      ))
    },.parallel = F
  )


}
