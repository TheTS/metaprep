
#' es_within
#'
#' @description Calculate effect size (md, d, g) from a within-subjects design
#'  or matched  groups (e.g. a crossover study).
#'
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param sd_diff .
#' @param r .
#' @param t .
#' @param p .
#' @param n Number of subjects
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#' @importFrom stats qt
#'
es_within <- function(x1, x2, sd1, sd2, sd_diff = NULL, r = NULL, t = NULL, p = NULL, n) {

  md <- x1 - x2

  if (!is.null(r)) {
    sd_diff <- sqrt((sd1 ^ 2 + sd2 ^ 2) - (2 * r * sd1 * sd2))
    se_md <- sd_diff / sqrt(n)
  } else if (!is.null(p)) {
    if (is.null(t)) {
      t <- qt(p / 2, df = n - 1, lower.tail = F)
    }

    se_md <- md / t
    sd_diff <- se_md * sqrt(n)

  } else
    stop("Missing info - p must be specified if r is missing.")

  variance_md <- (sd_diff ^ 2) / n

  if (is.null(r)) {
    r <- (sd1 ^ 2 + sd2 ^ 2 - sd_diff ^ 2) / (2 * sd1 * sd2)
  }

  sd_within <- sd_diff / sqrt(2 * (1 - r))
  d <- md / sd_within
  variance_d <- (1 / n) + (d ^ 2) / (2 * n) * (2 * (1 - r))
  se_d <- sqrt(variance_d)
  j <- 1 - (3 / (4 * (n - 1) - 1))
  g <- j * d
  variance_g <- j ^ 2 * variance_d
  se_g <- sqrt(variance_g)


  data.frame(md, variance_md, se_md, d, variance_d, se_d, g, variance_g, se_g)
}


#' es_between
#'
#' @description Calculate effect size (md, d, g) from a between-subjects design
#'  of independent groups (e.g. a parallel randomised controlled trial).
#'
#' @param x1 Mean of group 1
#' @param x2 Mean of group 2
#' @param sd1 Standard deviation of group 1
#' @param sd2 Standard deviation of group 2
#' @param n1 Number of subjects in group 1
#' @param n2 Number of subjects in group 2
#'
#' @return
#'
#' @export
#'
#' @examples
#'
es_between <- function(x1, x2, sd1, sd2, n1, n2) {

  md <- x1 - x2
  sPooled <- sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2) / (n1 + n2 - 2))
  variance_md <- ((n1 + n2) / (n1 * n2)) * sPooled ^ 2
  #var_diff_nonorm <- (sd1 ^ 2 / n1) + (sd2 ^ 2 / n2)
  se_md <- sqrt(variance_md)
  #se_diff_nonorm <- sqrt(var_diff_nonorm)
  d <- (x1 - x2) / sPooled
  variance_d <- ((n1 + n2) / (n1 * n2)) + d ^ 2 / (2 * (n1 + n2))
  se_d <- sqrt(variance_d)
  j <- 1 - (3 / (4 * ((n1 + n2) - 2) - 1))
  g <- j * d
  variance_g = j ^ 2 * variance_d
  se_g = sqrt(variance_g)

  data.frame(md, variance_md, se_md, d, variance_d, se_d, g, variance_g, se_g)
}



#' combine_within
#'
#' @description Combine effect sizes from the same within-subjects study. An
#' example is a crossover study with multiple intervention conditions.
#'
#' @param es A vector of effect sizes
#' @param var_es A vector of effect size variances
#' @param r The correlation
#'
#' @return
#'
#' @export
#'
#' @examples
#'
combine_within <- function (es, var_es, r) {

    n_studies <- length(es)
    mean_effect <- sum(es) / n_studies
    variance <- sum(var_es) / n_studies
    se <-  sum(sqrt(var_es)) / n_studies
    mean_variance <- (variance / n_studies) * (1 - r) + se ^ 2 * r

    res <- data.frame(mean_effect, mean_variance, se)
}



#' combine_between
#'
#' @description Combining mean and standard devation from independant subgroups
#' within a study. This is computing a composite score by using the summary data
#' from the subgroups to recreate the data for the study as a whole. This summary
#' data can then be later used to compute the effect size and variance.
#'
#' @param x A vector of means
#' @param sd A vector of standard deviations
#' @param n A vector of number of subjects
#'
#' @return
#'
#' @export
#'
#' @examples
#'
combine_between <- function(x, sd, n) {

  n_studies = (length(x))

  n1 = n[1]
  n2 = n[2]
  sd1 = sd[1]
  sd2 = sd[2]
  x1 = x[1]
  x2 = x[2]

  combined_n = n1 + n2
  combined_x = (n1 * x1 + n2 * x2) / (n1 + n2)
  combined_sd = sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2 + ((n1 * n2) / (n1 + n2)) * ((x1 - x2) ^ 2)) / (n1 + n2 - 1))
  #combined_sd = sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2 + (((n1 * n2)/(n1 + n2)) * (x1^2 + x2^2 - 2 * x1 * x2)))/(n1 + n2 - 1))

  res <- data.frame(combined_n, combined_x, combined_sd)

  if (n_studies > 2) {

    for (i in 3:n_studies) {

      n1 = res$combined_n
      n2 = n[i]
      sd1 = res$combined_sd
      sd2 = sd[i]
      x1 = res$combined_x
      x2 = x[i]

      combined_n = n1 + n2
      combined_x = (n1 * x1 + n2 * x2) / (n1 + n2)
      combined_sd = sqrt(((n1 - 1) * sd1 ^ 2 + (n2 - 1) * sd2 ^ 2 + ((n1 * n2) / (n1 + n2)) * ((x1 - x2) ^ 2)) / (n1 + n2 - 1))

      res <- data.frame(combined_n, combined_x, combined_sd)
    }

  }

  res
}





#' pool_between
#'
#' @description Takes a dataframe of means and standard deviations and returns
#' a dataframe with pooled means and standard deviations for between-subject
#' studies. This function is normally run before calculating effect sizes.
#'
#' @param data A dataframe
#' @param n.e Column name: number of subjects in experimental group
#' @param mean.e Column name: mean for experimental group
#' @param sd.e Column name: Standard deviation for experimental group
#' @param n.c Column name: number of subjects in control group
#' @param mean.c Column name: mean for control group
#' @param sd.c Column name: Standard deviation for control group
#' @param type Column name: type of study.
#' @param studlab Column name: study label
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#' @import dplyr
#' @importFrom rlang :=
#'
pool_between <- function(data, n.e, mean.e, sd.e, n.c, mean.c, sd.c, type, studlab) {

  n.e <- enquo(n.e)
  mean.e <- enquo(mean.e)
  sd.e <- enquo(sd.e)
  n.c <- enquo(n.c)
  mean.c <- enquo(mean.c)
  sd.c <- enquo(sd.c)
  type <- enquo(type)
  studlab <- enquo(studlab)

  data %>%
    group_by(!!studlab) %>%
     group_map(function(df, group) {

      if (nrow(df) > 1 & pull(df, !!type)[1] == 'between') {

        cb <- combine_between(x = pull(df, !!mean.e),
                              sd = pull(df, !!sd.e),
                              n = pull(df, !!n.e))


        res <- df[1, ] %>%
          mutate(!!quo_name(n.c) := pull(df, !!n.c)[1],
                 !!quo_name(n.e) := cb$combined_n,
                 !!quo_name(mean.c) := pull(df, !!mean.c)[1],
                 !!quo_name(sd.c) := pull(df, !!sd.c)[1],
                 !!quo_name(mean.e) := cb$combined_x,
                 !!quo_name(sd.e) :=  cb$combined_sd)

        return(cbind(group, res))
      }

      cbind(group, df)
    }) %>%
    bind_rows()
}




#' pool_within
#'
#' @description Pool...
#'
#' @param data .
#' @param md .
#' @param variance_md .
#' @param se_md .
#' @param d .
#' @param variance_d .
#' @param se_d .
#' @param g .
#' @param variance_g .
#' @param se_g .
#' @param type .
#' @param studlab .
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#' @importFrom rlang quo_is_null
#' @importFrom rlang :=
#' @import dplyr
#'
pool_within <- function(data,
                        md=NULL, variance_md=NULL, se_md=NULL,
                        d=NULL, variance_d=NULL, se_d=NULL,
                        g=NULL, variance_g=NULL, se_g=NULL, type, studlab) {

  md <- enquo(md)
  variance_md <- enquo(variance_md)
  se_md <- enquo(se_md)
  d <- enquo(d)
  variance_d <- enquo(variance_d)
  se_d <- enquo(se_d)
  g <- enquo(g)
  variance_g <- enquo(variance_g)
  se_g <- enquo(se_g)
  type <- enquo(type)
  studlab <- enquo(studlab)


  data %>%
    group_by(!!studlab) %>%
    group_map(function(df, group) {

      if (nrow(df) > 1 & pull(df, !!type)[1] == 'within') {

        res <- df[1, ] %>%
          mutate_at(vars(-!!quo_name(type)), ~NA)


        if (!quo_is_null(d)) {
          cw <- combine_within(es = pull(df, !!d),
                               var_es = pull(df, !!variance_d),
                               r = 0.5)

          res <- res %>%
            mutate(!!quo_name(d) := cw$mean_effect,
                   !!quo_name(variance_d) := cw$mean_variance,
                   !!quo_name(se_d) := cw$se)
        }

        if (!quo_is_null(g)) {
          cw <- combine_within(es = pull(df, !!g),
                               var_es = pull(df, !!variance_g),
                               r = 0.5)

          res <- res %>%
            mutate(!!quo_name(g) := cw$mean_effect,
                   !!quo_name(variance_g) := cw$mean_variance,
                   !!quo_name(se_g) := cw$se)
        }

        if (!quo_is_null(md)) {
          cw <- combine_within(es = pull(df, !!md),
                               var_es = pull(df, !!variance_md),
                               r = 0.5)

          res <- res %>%
            mutate(!!quo_name(md) := cw$mean_effect,
                   !!quo_name(variance_md) := cw$mean_variance,
                   !!quo_name(se_md) := cw$se)
        }

        return(cbind(group, res))

      }

      cbind(group, df)
    }) %>%
    bind_rows()
}



#' calculate_es
#'
#' @description Calculate...
#'
#' @param data .
#' @param n.e .
#' @param mean.e .
#' @param sd.e .
#' @param n.c .
#' @param mean.c .
#' @param sd.c .
#' @param type .
#' @param studlab .
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#' @import dplyr
#' @importFrom rlang :=
#'
calculate_es <- function(data, n.e, mean.e, sd.e, n.c, mean.c, sd.c, type, studlab) {

  n.e <- enquo(n.e)
  mean.e <- enquo(mean.e)
  sd.e <- enquo(sd.e)
  n.c <- enquo(n.c)
  mean.c <- enquo(mean.c)
  sd.c <- enquo(sd.c)
  type <- enquo(type)
  studlab <- enquo(studlab)

  result <- data.frame()

  for (i in 1:nrow(data)) {
    df1 <- data[i, ]

    if (pull(df1, !!type)[1] == 'between') {

      x <- es_between(x1 = pull(df1, !!mean.e),
                      x2 = pull(df1, !!mean.c),
                      sd1 = pull(df1, !!sd.e),
                      sd2 = pull(df1, !!sd.c),
                      n1 = pull(df1, !!n.e),
                      n2 = pull(df1, !!n.c)) %>%
        mutate(!!quo_name(type) := pull(df1, !!type))

      result <- rbind(result, x)

    } else {
      x <- es_within(
        x1 = pull(df1, !!mean.e),
        x2 = pull(df1, !!mean.c),
        sd1 = pull(df1, !!sd.e),
        sd2 = pull(df1, !!sd.c),
        n = pull(df1, !!n.e),
        r = 0.5) %>%
        mutate(!!quo_name(type) := pull(df1, !!type))

      result <- rbind(result, x)
    }

  }

  cbind(study = pull(data, !!studlab), result)
}


#' pool_studies
#'
#' @description Takes a dataframe of means, standard deviations, number of
#' subjects in each group, and study type (either 'within' or 'between') and
#' returns a dataframe of study effect sizes (that are pooled across studies
#' where necessary). This can then be used for meta analysis.
#'
#' @param data .
#' @param n.e .
#' @param mean.e .
#' @param sd.e .
#' @param n.c .
#' @param mean.c .
#' @param sd.c .
#' @param type .
#' @param studlab .
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#' @import dplyr
#'
pool_studies <- function(data, n.e, mean.e, sd.e, n.c, mean.c, sd.c, type, studlab) {

  n.e <- enquo(n.e)
  mean.e <- enquo(mean.e)
  sd.e <- enquo(sd.e)
  n.c <- enquo(n.c)
  mean.c <- enquo(mean.c)
  sd.c <- enquo(sd.c)
  type <- enquo(type)
  studlab <- enquo(studlab)


  data %>%
    pool_between(n.e = !!n.e,
                 mean.e = !!mean.e,
                 sd.e = !!sd.e,
                 n.c = !!n.c,
                 mean.c = !!mean.c,
                 sd.c = !!sd.c,
                 type = !!type,
                 studlab = !!studlab) %>%

    calculate_es(n.e = !!n.e,
                 mean.e = !!mean.e,
                 sd.e = !!sd.e,
                 n.c = !!n.c,
                 mean.c = !!mean.c,
                 sd.c = !!sd.c,
                 type = !!type,
                 studlab = !!studlab) %>%

    pool_within(md = .data$md,
                variance_md = .data$variance_md,
                se_md = .data$se_md,
                d = .data$d,
                variance_d = .data$variance_d,
                se_d = .data$se_d,
                g = .data$g,
                variance_g = .data$variance_g,
                se_g = .data$se_g,
                type = !!type,
                studlab = !!studlab)
}


