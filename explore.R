library(tidyverse)
library(scales)

years = 2010:2019
regions = 1:10


data = (
  cross(list(year=years, region=regions))
  %>% map_chr(~ file.path(
                "data", "ILINet",
                .x $ year, paste0(.x $ region, ".csv")
              ))
  %>% keep(file.exists)
  %>% map_df(read_csv, col_types=cols(), id="fname")
  %>% mutate(
    region=basename(fname) %>% str_sub(1, -4) %>% as.integer,
    year=dirname(fname) %>% basename %>% as.integer
  )
  %>% select(-fname)
)


# Plot how well trajectories fit
(trajs
  %>% group_by(year, region, t)
  %>% summarise(I=mean(I))
  %>% rename(week=t)
  %>% inner_join(data, by=c("year", "region", "week"))
  %>% select(week, dS, I, year, region)
  %>% gather(var, val, -year, -region, -week)
  %>% ggplot(aes(week, val, color=var))
  + facet_grid(year ~ region, scales='free', label=label_both)
  + theme_minimal()
  + theme(legend.position="bottom")
  + geom_line()
)

## (data
##   %>% filter(year==2017)
##   %>% ggplot(aes(week, dS, color=factor(region)))
##   + geom_line()
##   + theme_minimal()
##   + theme(legend.position="bottom")
## )

peaks = (
  data
  %>% group_by(year, region)
  %>% arrange(week)
  %>% slice(which.max(dS))
  %>% rename(peak_t=week, peak_dS=dS)
  %>% select(-Nmax)
)

results = (
  cross(list(year=years, region=regions))
  %>% map_chr(~ file.path(
                "results", "ILINet", "estimates",
                .x $ year, paste0(.x $ region, ".csv")
              ))
  %>% keep(file.exists)
  %>% map_df(read_csv, col_types=cols(), id="fname")
  %>% mutate(
    region=basename(fname) %>% str_sub(1, -4) %>% as.integer,
    year=dirname(fname) %>% basename %>% as.integer
  )
  %>% select(-fname)
)

fisher_ub <- function(N, beta, gamma, x)
  x^2 / ((N - x) * N^2 * (N - x + N * gamma / beta))

cum_fisher_ub <- function(N, beta, gamma, x)
  (1:x)^2 / ((N - 1:x) * N^2 * (N - 1:x + N * gamma / beta))

thm1_samples = function(N, beta, gamma, thresh=1)
  which.max(((1 / cumsum(fisher_ub(N, beta, gamma, 1:N))) / (N ^ 2)) < sqrt(thresh))

## fisher_ub <- function(N, beta, gamma, x)
##   9 * x^2 / (0.5 * (0.5 + gamma / beta) * N^4)


z = (1 / cumsum(fisher_ub(43817, 0.124, 0.24, 1:43817))) / (43817 ^ 2)

which.max(z < 1)

thm1_samples(43817, 0.124, 0.24)

final_results <- (
  results
  %>% group_by(year, region)
  %>% slice(which.max(t))
  %>% mutate(
        m_star=thm1_samples(N, beta, gamma),
        m_star_50=thm1_samples(N, beta, gamma, 0.5)
      )
  %>% rename(Ntrue=N)
  %>% select(
        region, year, m_star,
        m_star_50, Ntrue, beta, gamma)
)

to_plot = (final_results
  %>% inner_join(peaks, by=c("year", "region"))
  %>% inner_join(data, by=c("year", "region"))
  %>% group_by(year, region, peak_t)
  %>% arrange(week)
  %>% mutate(m=cumsum(dS))
  %>% summarise(
        const_frac_t=which.max(m - min(m) > Ntrue / 4),
        mstar_t=which.max(m - min(m) > m_star )
        )
)

(to_plot
  %>% gather(var, val, -year, -region)
  %>% ggplot(aes(val, color=var))
  + stat_ecdf()
  + scale_y_continuous(label=percent, name="Cum %")
  + theme_minimal()
  + facet_wrap(~ year, label=label_both)
  + theme(legend.position="bottom")
  + coord_cartesian(xlim=c(0, 35))
)

to_plot %>% ungroup %>% summarise(median(mstar_t / const_frac_t))



# Read trajectories
trajs = (
  cross(list(year=years, region=regions))
  %>% map_chr(~ file.path(
                "results", "ILINet", "trajectories",
                .x $ year, paste0(.x $ region, ".csv")
              ))
  %>% keep(file.exists)
  %>% map_df(read_csv, col_types=cols(), id="fname")
  %>% mutate(
    region=basename(fname) %>% str_sub(1, -4) %>% as.integer,
    year=dirname(fname) %>% basename %>% as.integer
  )
  %>% select(-fname)
)

# Add dS
trajs_with_dS = (
  trajs
  %>% group_by(region, year, seed)
  %>% arrange(t)
  %>% mutate(dS=I - lag(I))
)


fisher_ub <- function(N, beta, gamma, x)
  x^2 / ((N - x) * N^2 * (N - x + N * gamma / beta))


data.frame(
  Ntrue =
)

trajs_agg = (
  trajs_with_dS
  %>% inner_join(final_results, by=c("region", "year"))
  %>% mutate(C=I + R)
  %>% mutate(fisher=C^2 / (Ntrue - C) / Ntrue^2 / (Ntrue - C + Ntrue * gamma / beta)) %>% mutate(fisher=map(C, ~ cum_fisher_ub(Ntrue, beta, gamma, .x)))
  %>% group_by(year, region, t, Ntrue)
  %>% summarise(dS=mean(dS), m_star_t=mean(fisher))
)

# Aggregate over seeds

trajs_agg = (trajs_with_dS
  %>% inner_join(final_results, by=c("region", "year"))
  %>% mutate(C=I + R)
  %>% group_by(year, region, Ntrue, m_star, m_star_50, seed)
  %>% summarise(
        m_star_t=which.max(C > m_star),
        m_star_50_t=which.max(C > m_star_50),
        m_peak_t=which.max(I),
        m_peak=C[which.max(I)],
        m_tot=sum(I)
      )
  %>% group_by(year, region, Ntrue, m_star, m_star_50)
  %>% summarise(
        m_star_t=mean(m_star_t),
        m_star_50_t=mean(m_star_50_t),
        m_peak_t=mean(m_peak_t),
        m_peak=mean(m_peak),
        m_tot=mean(m_tot)
      )
)







  ## %>% mutate(fisher=C^2 / (Ntrue - C) / Ntrue^2 / (Ntrue - C + Ntrue * gamma / beta))
  ## %>% mutate(fisher=map(C, ~ cum_fisher_ub(Ntrue, beta, gamma, .x)))
  %>% group_by(year, region, t, Ntrue)
  %>% summarise(dS=mean(dS), fisher=mean(fisher))
  %>% ungroup
  %>% group_by(year, region, Ntrue)
  %>% arrange(t)
  ## %>% mutate(rel_err=1 / cumsum(fisher) / max(Ntrue)^2)
)

m_stars = (
  trajs_agg
  %>% group_by(region, year, Ntrue)
  %>% arrange(t)
  %>% mutate(m=cumsum(coalesce(dS, 0)))
  %>% mutate(fisher=map(C, ~ cum_fisher_ub(Ntrue, beta, gamma, .x)))
  %>% summarise(
        m_star_t=which.max(rel_err < 1),
        m_star=m[which.max(rel_err < 1)],
        m_peak_t=which.max(dS),
        m_peak=m[which.max(dS)],
      )
)

(m_stars %>% ggplot(aes(m_star_t, m_peak_t))
  + geom_point()
  + theme_minimal()
  + theme(legend.position="bottom")
  + geom_abline(linetype='dashed', alpha=0.4)
)


library(latex2exp)
library(viridis)

leftovers = (trajs
  %>% group_by(region, year, seed)
  %>% filter(t > 52)
  %>% summarise(n=sum(I))
  %>% group_by(region, year)
  %>% summarise(m_leftover=mean(n))
  %>% mutate(year=year + 1)
)

toplot = (data
  %>% filter(year > 2010)
  %>% inner_join(trajs_agg, by=c("year", "region"))
  %>% inner_join(leftovers, by=c("year", "region"))
  ## %>% inner_join(m_stars, by=c("year", "region"))
  ## %>% inner_join(final_results %>% select(region, year, m_thm1=m_star)
  ##              , by=c("year", "region"))
  %>% group_by(year, region)
  %>% arrange(week)
  %>% mutate(m=cumsum(dS) - m_leftover)
  %>% summarise(
        m_star_t=which.max(m > m_star),
        m_star_50_t=which.max(m > m_star_50),
        m_peak_t=coalesce(which.max(m > m_peak / 4), 52),
        peak_t=which.max(dS),
      )
)

(toplot
  %>% ggplot(aes(m_star_t / peak_t))
  + theme_minimal()
  + theme(legend.position="bottom")
  + stat_ecdf()
  + scale_y_continuous(label=percent, name="Cum % Region-Years")
  + scale_x_continuous(label=percent, name="Weeks to RelError = 1 / Weeks to Peak")
)


(toplot
  %>% ggplot(aes(m_star_t, peak_t))
  + theme_minimal()
  + theme(legend.position="bottom")
  + geom_point()
  + geom_abline(linetype='dashed', alpha=0.4)
)


(toplot
  %>% ggplot(aes(m_thm1_t))
  + stat_ecdf()
  + scale_y_continuous(label=percent, name="Cum %")
  + theme_minimal()
  + theme(legend.position="bottom")
)


library(patchwork)

rel = (toplot
  %>% ggplot(aes(m_star_t, y=..count.. / sum(..count..)))
  + geom_histogram(binwidth=1)
  + scale_y_continuous(label=percent, name="% Instances")
  + theme_minimal()
  + theme(legend.position="bottom")
  + scale_x_continuous(
      breaks=seq(0, 52, 2),
      name=("Weeks to RelError = 1")
    )
  + coord_cartesian(ylim=c(0, 0.25))
  )
pk = (toplot
  %>% ggplot(aes(m_peak_t - m_star_t, y=..count.. / sum(..count..)))
  + geom_histogram(binwidth=1)
  + scale_y_continuous(label=percent, name=NULL)
  + theme_minimal()
  + theme(legend.position="bottom")
  + scale_x_continuous(
      breaks=seq(-20, 52, 1),
      name=("Weeks from RelError = 1 to Peak")
    )
  + coord_cartesian(ylim=c(0, 0.25))
  )
rel | pk


(toplot
  ## %>% gather(var, val, -year, -region)
  ## %>% group_by(var, val)
  ## %>% summarise(n=n())
  ## %>% group_by(var)
  ## %>% mutate(pct=n / sum(n))
  %>% group_by(m_star_t, m_peak_t)
  %>% summarise(n=n())
  %>% ungroup
  %>% mutate(pct=round(n * 100 / sum(n)) / 100)
  %>% ggplot(aes(m_star_t, m_peak_t, fill=pct))
  + geom_tile()
  + geom_text(aes(label=scales::percent(pct)), size=3,
              color="gray", fontface="bold")
  + geom_abline(linetype='dashed', alpha=0.4)
  + scale_fill_viridis()
  + scale_x_continuous(
      breaks=seq(0, 52, 2),
      name=TeX("Weeks to RelError = $o(1)$")
    )
  + scale_y_continuous(
      breaks=seq(0, 52, 2),
      name="Weeks to Peak"
    )
  + theme_minimal()
  + guides(fill=FALSE)
  ## + scale_y_continuous(label=percent, name="Cum % Instances")
  ## + scale_y_continuous(label=percent, name="Cum %")
)

(
  toplot
  %>% select(region, year, m_star_t, m_peak_t)
  %>% gather(var, val, -region, -year)
  %>% mutate(var=case_when(
        var=="m_peak_t" ~ "Weeks to Ctot / 8 infections",
        var=="m_star_t" ~ "Weeks to RelError bound = 1",
        TRUE ~ var
        )
      )
  %>% ggplot(aes(val, linetype=var))
  + stat_ecdf()
  + scale_y_continuous(label=percent, name="Cum % Instances")
  + xlab("Weeks")
  + theme_minimal()
  + theme(legend.position="bottom", legend.title=element_blank())
)



(data
  %>% filter(year > 2010)
  %>% inner_join(m_stars, by=c("year", "region"))
  %>% group_by(year, region)
  %>% arrange(week)
  %>% mutate(m=cumsum(dS))
  %>% summarise(
        m_star_t=which.max(m > m_star),
        m_peak_t=coalesce(which.max(m > m_peak), 52),
        )
  ## %>% gather(var, val, -year, -region)
  ## %>% group_by(var, val)
  ## %>% summarise(n=n())
  ## %>% group_by(var)
  ## %>% mutate(pct=n / sum(n))
  %>% group_by(m_star_t, m_peak_t)
  %>% summarise(n=n())
  %>% ungroup
  %>% mutate(pct=round(n * 100 / sum(n)) / 100)
  %>% ggplot(aes(m_star_t, m_peak_t, size=pct))
  + geom_point()
  ## + geom_text(aes(label=scales::percent(pct)), size=3,
  ##             color="gray", fontface="bold")
  + geom_abline(linetype='dashed', alpha=0.4)
  + scale_fill_viridis()
  + scale_x_continuous(
      breaks=seq(0, 52, 2),
      name=TeX("Weeks to RelError = $o(1)$")
    )
  + scale_y_continuous(
      breaks=seq(0, 52, 2),
      name="Weeks to Peak"
    )
  + theme_minimal()
  + guides(fill=FALSE)
  ## + scale_y_continuous(label=percent, name="Cum % Instances")
  ## + scale_y_continuous(label=percent, name="Cum %")
)


(final_results
  %>% ggplot(aes(Ntrue))
  + theme_minimal()
  + theme(legend.position="bottom")
  + stat_ecdf()
  + scale_y_continuous(label=percent, name="Cum % Region-Years")
  + scale_x_log10(name="Estimated N")
)



# Amazon
x = (
  list.files('results/Amazon/estimates', full.names=TRUE)
  %>% map_df(read_csv, col_types=cols(), id="fname")
)

(x
  %>% filter(N < 1e5 -10)
  %>% ggplot(aes(N))
+ stat_ecdf()
+ scale_y_continuous(label=percent, name="Cum % Products")
+ theme_minimal()
+ theme(legend.position="bottom")
+ scale_x_log10(name="Estimated N")
)

x %>% filter(N < 1e5 -10) %>% summarise(mean(N))
