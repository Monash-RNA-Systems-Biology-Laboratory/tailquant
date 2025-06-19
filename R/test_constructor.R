
make_design <- function(samples, groups, batches=NULL) {
    assertthat::assert_that(is.factor(groups))
    assertthat::assert_that(is.null(batches) || is.factor(batches))
    
    n <- length(samples)
    if (!is.null(batches) && length(unique(batches)) > 1) {
        batch_design <- contr.treatment(levels(batches))[as.numeric(batches),,drop=FALSE]
    } else {
        batch_design <- matrix(0, nrow=n, ncol=0)
    }
    
    group_design <- contr.treatment(levels(groups), contrasts=FALSE)[as.numeric(groups),,drop=FALSE]
    design <- cbind(group_design, batch_design)
    rownames(design) <- samples
    
    design
}

#' @export
make_test <- function(samples, groups, group1, group2, batches=NULL, name=NULL, title=NULL) {
    assertthat::assert_that(length(group1) == length(group2))
    
    design <- make_design(samples, groups, batches)
    coefs <- colnames(design)
    
    nc <- length(group1)
    contrasts <- matrix(0, nrow=length(coefs), ncol=nc)
    for(i in seq_len(nc)) {
        contrasts[ coefs == group1[i], i ] <- -1
        contrasts[ coefs == group2[i], i ] <- 1
    }
    colnames(contrasts) <- paste0(group2,"-",group1)
    
    test <- list(
        title = title %||% paste0(group1, " to ", group2),
        design = design,
        contrasts = contrasts)
    
    result <- list(test)
    names(result) <- name %||% paste0(group1,"_to_",group2, collapse=" or ")
    result
}

make_test_oneway_anova <- function(samples, groups, batches=NULL) {
    assertthat::assert_that(is.factor(groups))
    assertthat::assert_that(is.null(batches) || is.factor(batches))
    
    levels <- levels(groups)
    n <- length(levels)
    
    make_test(
        samples, groups, rep(levels[1],n-1), levels[-1], batches, 
        name="any", title="Any change")
}

#' @export
make_tests_oneway <- function(samples, groups, batches=NULL) {
    assertthat::assert_that(is.factor(groups))
    assertthat::assert_that(is.null(batches) || is.factor(batches))
    
    levels <- levels(groups)
    n <- length(levels)
    
    result <- list()
    if (n > 2) {
        result <- c(result, 
            make_test_oneway_anova(samples, groups, batches))
    }
    
    for(i in seq(1,n-1)) {
        for(j in seq(i+1,n)) {
            result <- c(result, make_test(samples, groups, levels[i], levels[j], batches))
        }
    }
    
    result
}

make_tests_twoway_helper <- function(samples, groups1, groups2, batches=NULL) {
    assertthat::assert_that(is.factor(groups1))
    assertthat::assert_that(is.factor(groups2))
    assertthat::assert_that(is.null(batches) || is.factor(batches))
    
    result <- list()
    
    groups <- forcats::fct_cross(groups1, groups2, sep="_")
    
    levels <- levels(groups1)
    n <- length(levels)
    
    for(level2 in levels(groups2)) {
        this_levels <- paste0(levels,"_",level2)
        if (n > 2) {
            result <- c(result,
                make_test(samples, groups,
                    rep(this_levels[1], n-1), this_levels[-1], batches,
                    name=paste0("any_",level2), title=paste0("Any change within ", level2)))
        }
        
        for(i in seq(1,n-1)) {
            for(j in seq(i+1,n)) {
                result <- c(result, 
                    make_test(samples, groups, 
                        this_levels[i], 
                        this_levels[j], 
                        batches))
            }
        }
    }
    
    result
}

make_test_twoway_interaction <- function(samples, groups1, groups2, batches=NULL) {
    assertthat::assert_that(is.factor(groups1))
    assertthat::assert_that(is.factor(groups2))
    assertthat::assert_that(is.null(batches) || is.factor(batches))
    
    groups <- forcats::fct_cross(groups1, groups2, sep="_")
    levels <- levels(groups)
    
    levels1 <- levels(groups1)
    n1 <- length(levels1)
    levels2 <- levels(groups2)
    n2 <- length(levels2)
    
    design <- make_design(samples, groups, batches)
    coefs <- colnames(design)
    
    contrasts <- matrix(0, nrow=length(coefs), ncol=(n1-1)*(n2-1))
    colnames(contrasts) <- rep("", ncol(contrasts))
    
    j <- 1
    for(i1 in seq(2,n1)) {
        for(i2 in seq(2,n2)) {
            contrasts[coefs == paste0(levels1[ 1],"_",levels2[ 1]), j] <-  1
            contrasts[coefs == paste0(levels1[i1],"_",levels2[ 1]), j] <- -1
            contrasts[coefs == paste0(levels1[ 1],"_",levels2[i2]), j] <- -1
            contrasts[coefs == paste0(levels1[i1],"_",levels2[i2]), j] <-  1
            colnames(contrasts)[j] <- paste0("int_",levels1[i1],"_",levels2[i2])
            j <- j+1
        }
    }
    
    test <- list(
        title = "Any interaction",
        design = design,
        contrasts = contrasts)
    
    result <- list(test)
    names(result) <- "any_interaction"
    result
}

#' @export
make_tests_twoway <- function(samples, groups1, groups2, batches=NULL) {
    c(
        make_test_oneway_anova(samples, forcats::fct_cross(groups1,groups2,sep="_"), batches),
        make_test_twoway_interaction(samples, groups1, groups2, batches),
        make_tests_twoway_helper(samples, groups1, groups2, batches),
        make_tests_twoway_helper(samples, groups2, groups1, batches))
}