
theme <- function() {
    ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank())
}

#' @export
plot_tail_detail <- function(reads, min_tail=0) {
    # Accept data in many forms
    if (!identical(class(reads), "list")) {
        reads <- list(reads)
    }
    
    df <- reads |>
        purrr::map(dplyr::collect) |>
        dplyr::bind_rows() |>
        dplyr::count(tail_start, tail) |>
        dplyr::filter(tail >= .env$min_tail)
    
    ggplot2::ggplot(df)+
        ggplot2::aes(x=tail_start,y=tail,fill=sqrt(n))+
        ggplot2::geom_hline(yintercept=min_tail-0.5)+
        ggplot2::geom_tile(width=1, height=1)+
        ggplot2::scale_fill_viridis_c(labels=\(i) scales::comma(i*i))+
        ggplot2::coord_fixed(expand=FALSE, xlim=c(0,max(df$tail_start)+1), ylim=c(0,max(df$tail)+1))+
        ggplot2::labs(x="Tail start",y="Tail length",fill="Count")+
        theme() +
        ggplot2::theme(panel.background = ggplot2::element_rect("#666666"))
}

#' @export
plot_km_survival <- function(kms, names="Data", min_tail, max_tail=NULL, cpm=FALSE, lib_sizes=NULL) {
    kms <- ensure_list(kms)
    
    if (is.null(max_tail))
         max_tail <- max(purrr::map_dbl(kms, \(item) max(min_tail, item$tail)))
    
    unit <- "Reverse cumulative distribution"
    scale <- rep(1, length(kms))
    if (cpm) {
        unit <- "Reverse cumulative CPM"
        for(i in seq_along(kms)) {
            scale[i] <- kms[[i]]$active_n[1] * 1e6 / lib_sizes[i]
            kms[[i]]$prop_before <- kms[[i]]$prop_before * scale[i]
            kms[[i]]$prop_after  <- kms[[i]]$prop_after  * scale[i]
        }
    }
    
    names <- forcats::fct_inorder(names)
    df <- dplyr::tibble(name=.env$names, km=.env$kms) |>
         tidyr::unnest("km") |>
         tidyr::pivot_longer(c(prop_before, prop_after), names_to="what", values_to="prop") |>
         dplyr::select(name, tail, prop)
    df <- dplyr::bind_rows(dplyr::tibble(name=.env$names, tail=min_tail,prop=scale), df)
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail,y=prop,color=name,group=name) + 
        ggplot2::geom_path() + 
        ggplot2::coord_cartesian(xlim=c(min_tail,max_tail), ylim=c(0,max(scale))) +
        ggplot2::labs(x="Tail length", y=unit, color="") + 
        theme()
}

#' @export
plot_km_density <- function(kms, names="Data", min_tail=0, max_tail=NULL, step=1, cpm=FALSE, lib_sizes=NULL) {
    kms <- ensure_list(kms)
    
    unit <- "Proportion"
    if (cpm) {
        unit <- "CPM"
        for(i in seq_along(kms))
            kms[[i]]$prop_died <- kms[[i]]$prop_died * kms[[i]]$active_n[1] * 1e6 / lib_sizes[i]
    }
    
    kms <- purrr::map(kms, km_complete, min_tail=min_tail, max_tail=max_tail)
    
    df <- dplyr::tibble(name=forcats::fct_inorder(.env$names), km=.env$kms) |>
        tidyr::unnest("km") |>
        dplyr::mutate(tail = (floor((tail-min_tail)/step)+0.5)*step+min_tail) |>
        dplyr::summarize(prop_died=sum(prop_died), .by=c(name,tail))
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail,y=prop_died,color=name,group=name) + 
        #geom_path(aes(y=n_event/sum(n_event)), color="#ff0000") +
        ggplot2::geom_path() +
        ggplot2::coord_cartesian(ylim=c(0,max(df$prop_died))) +
        ggplot2::labs(x="Tail length", y=unit, color="") +
        theme()
}


#' @export
plot_km_density_ridgeline <- function(kms, names="Data", min_tail=0, max_tail=NULL, step=1, cpm=FALSE, lib_sizes=NULL) {
    kms <- ensure_list(kms)
    
    unit <- "Proportion"
    if (cpm) {
        unit <- "CPM"
        for(i in seq_along(kms))
            kms[[i]]$prop_died <- kms[[i]]$prop_died * kms[[i]]$active_n[1] * 1e6 / lib_sizes[i]
    }
    
    kms <- purrr::map(kms, km_complete, min_tail=min_tail, max_tail=max_tail)
    
    df <- dplyr::tibble(name=forcats::fct_inorder(.env$names), km=.env$kms) |>
        tidyr::unnest("km") |>
        dplyr::mutate(tail = (floor((tail-min_tail)/step)+0.5)*step+min_tail) |>
        dplyr::summarize(prop_died=sum(prop_died), .by=c(name,tail))
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail,y=forcats::fct_rev(name),height=prop_died,fill=name) + 
        ggridges::geom_ridgeline(scale=1.5/max(df$prop_died), alpha=0.5) +
        ggplot2::coord_cartesian(ylim=c(1,length(names) + 1.5)) +
        ggplot2::labs(x="Tail length", y=unit, fill="") +
        ggplot2::guides(fill="none") +
        theme() +
        ggplot2::theme(panel.grid.major.y = ggplot2::element_line(size = 0.5))
}


#' @export
plot_km_density_heatmap <- function(kms, names="Data", min_tail=0, max_tail=NULL, step=1, normalize_max=TRUE, cpm=FALSE, lib_sizes=NULL) {
    kms <- ensure_list(kms)
    
    unit <- "Proportion"
    if (cpm && !normalize_max) {
        unit <- "CPM"
        for(i in seq_along(kms))
            kms[[i]]$prop_died <- kms[[i]]$prop_died * kms[[i]]$active_n[1] * 1e6 / lib_sizes[i]
    }
    
    if (is.null(max_tail))
        max_tail <- purrr::map_dbl(kms, \(km) max(km$tail)) |> max()
    
    kms <- purrr::map(kms, km_complete, min_tail=min_tail, max_tail=max_tail)
    
    df <- dplyr::tibble(name=forcats::fct_inorder(.env$names), km=.env$kms) |>
        tidyr::unnest("km") |>
        dplyr::mutate(tail = (floor((tail-min_tail)/step)+0.5)*step+min_tail) |>
        dplyr::summarize(prop_died=sum(prop_died), .by=c(name,tail))
    
    if (normalize_max) {
        unit <- "Proportion\n of maximum"
        df <- df |>
            dplyr::mutate(prop_died=prop_died / max(prop_died), .by=c(name))
    }
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail,y=forcats::fct_rev(name), fill=prop_died) + 
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::coord_cartesian(expand=FALSE) +
        ggplot2::labs(x="Tail length", y="", fill=unit) +
        theme()
}