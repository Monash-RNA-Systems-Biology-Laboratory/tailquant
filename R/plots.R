
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
        dplyr::count(genomic_bases, tail) |>
        dplyr::filter(tail >= .env$min_tail)
    
    ggplot2::ggplot(df)+
        ggplot2::aes(x=genomic_bases,y=tail,fill=sqrt(n))+
        ggplot2::geom_hline(yintercept=min_tail-0.5)+
        ggplot2::geom_tile(width=1, height=1)+
        ggplot2::scale_fill_viridis_c(labels=\(i) scales::comma(i*i))+
        ggplot2::coord_fixed(expand=FALSE, xlim=c(0,max(df$genomic_bases)+1), ylim=c(0,max(df$tail)+1))+
        ggplot2::labs(x="Genomic bases",y="Tail length",fill="Count")+
        theme() +
        ggplot2::theme(panel.background = ggplot2::element_rect("#666666"))
}

#' @export
plot_km_survival <- function(kms, names="Data", min_tail, max_tail=NULL) {
    kms <- ensure_list(kms)
    
    if (is.null(max_tail))
         max_tail <- max(purrr::map_dbl(kms, \(item) max(min_tail, item$tail)))
    
    names <- forcats::fct_inorder(names)
    df <- dplyr::tibble(name=.env$names, km=.env$kms) |>
         tidyr::unnest("km") |>
         tidyr::pivot_longer(c(prop_before, prop_after), names_to="what", values_to="prop") |>
         dplyr::select(name, tail, prop)
    df <- dplyr::bind_rows(dplyr::tibble(name=.env$names, tail=min_tail,prop=1), df)
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail,y=prop,color=name,group=name) + 
        ggplot2::geom_path() + 
        ggplot2::coord_cartesian(xlim=c(min_tail,max_tail), ylim=c(0,1)) +
        ggplot2::labs(x="Tail length", y="Reverse cumulative distribution") + 
        theme()
}

#' @export
plot_km_density <- function(kms, names="Data", min_tail=0, max_tail=NULL, step=1) {
    kms <- ensure_list(kms) |>
        purrr::map(km_complete, min_tail=min_tail, max_tail=max_tail)
    
    df <- dplyr::tibble(name=forcats::fct_inorder(.env$names), km=.env$kms) |>
        tidyr::unnest("km") |>
        dplyr::mutate(tail = (floor((tail-min_tail)/step)+0.5)*step+min_tail) |>
        dplyr::summarize(prop_died=sum(prop_died), .by=c(name,tail))
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail,y=prop_died,color=name,group=name) + 
        #geom_path(aes(y=n_event/sum(n_event)), color="#ff0000") +
        ggplot2::geom_path() +
        ggplot2::coord_cartesian(ylim=c(0,max(df$prop_died))) +
        ggplot2::labs(x="Tail length", y="Proportion") +
        theme()
}


#' @export
plot_km_density_heatmap <- function(kms, names="Data", min_tail=0, max_tail=NULL, step=1, normalize_max=TRUE) {
    if (is.null(max_tail))
        max_tail <- purrr::map_dbl(kms, \(km) max(km$tail)) |> max()
    
    kms <- ensure_list(kms) |>
        purrr::map(km_complete, min_tail=min_tail, max_tail=max_tail)
    
    df <- dplyr::tibble(name=forcats::fct_inorder(.env$names), km=.env$kms) |>
        tidyr::unnest("km") |>
        dplyr::mutate(tail = (floor((tail-min_tail)/step)+0.5)*step+min_tail) |>
        dplyr::summarize(prop_died=sum(prop_died), .by=c(name,tail))
    
    legend_title <- "Proportion"
    if (normalize_max) {
        legend_title <- "Proportion\n of maximum"
        df <- df |>
            dplyr::mutate(prop_died=prop_died / max(prop_died), .by=c(name))
    }
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail,y=forcats::fct_rev(name), fill=prop_died) + 
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::coord_cartesian(expand=FALSE) +
        ggplot2::labs(x="Tail length", y="", fill=legend_title) +
        theme()
}