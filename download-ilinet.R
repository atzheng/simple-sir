library(tidyverse)

#' Retrieve ILINet Surveillance Data
#' From https://github.com/hrbrmstr/cdcfluview/blob/60ffe50553371148e962daf5559fce8270ce3fdd/R/ilinet.r#L34
#'
#' The CDC FluView Portal provides in-season and past seasons' national, regional,
#' and state-level outpatient illness and viral surveillance data from both
#' ILINet (Influenza-like Illness Surveillance Network) and WHO/NREVSS
#' (National Respiratory and Enteric Virus Surveillance System).
#'
#' This function retrieves current and historical ILINet surveillance data for
#' the identified region.
#'
#' @md
#' @param region one of "`national`", "`hhs`", "`census`", or "`state`"
#' @param years a vector of years to retrieve data for (i.e. `2014` for CDC
#'        flu season 2014-2015). CDC has data for this API going back to 1997.
#'        Default value (`NULL`) means retrieve **all** years. NOTE: if you
#'        happen to specify a 2-digit season value (i.e. `57` == 2017-2018)
#'        the function is smart enough to retrieve by season ID vs convert that
#'        to a year.
#' @references
#' - [CDC FluView Portal](https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html)
#' - [ILINet Portal](https://wwwn.cdc.gov/ilinet/) (Login required)
#' - [WHO/NREVSS](https://www.cdc.gov/surveillance/nrevss/index.html)
#' @export
#' @examples
#' national_ili <- ilinet("national", years = 2017)
#' \dontrun{
#' hhs_ili <- ilinet("hhs")
#' census_ili <- ilinet("census")
#' state_ili <- ilinet("state")
#'
#' all_ili <- suppressWarnings(
#'   suppressMessages(purrr::map_df(c("national", "hhs", "census", "state"), ilinet)))
#' }
ilinet <- function(region = c("national", "hhs", "census", "state"), years = NULL) {

  #region="national"; years=1997:2018

  region <- match.arg(tolower(region), c("national", "hhs", "census", "state"))

  meta <- jsonlite::fromJSON("https://gis.cdc.gov/grasp/flu2/GetPhase02InitApp?appVersion=Public")

  list(
    AppVersion = "Public",
    DatasourceDT = list(list(ID = 1, Name = "ILINet")),
    RegionTypeId = .region_map[region]
  ) -> params

  params$SubRegionsDT <- switch(region,
    national = {
      list(list(ID = 0, Name = ""))
    },
    hhs = {
      lapply(1:10, function(i) list(ID = i, Name = as.character(i)))
    },
    census = {
      lapply(1:9, function(i) list(ID = i, Name = as.character(i)))
    },
    state = {
      lapply(1:59, function(i) list(ID = i, Name = as.character(i)))
    }
  )

  available_seasons <- sort(meta$seasons$seasonid)

  if (is.null(years)) { # ALL YEARS
    years <- available_seasons
  } else { # specified years or seasons or a mix

    years <- as.numeric(years)
    years <- ifelse(years > 1996, years - 1960, years)
    years <- sort(unique(years))
    years <- years[years %in% available_seasons]

    if (length(years) == 0) {
      years <- rev(sort(meta$seasons$seasonid))[1]
      curr_season_descr <- meta$seasons[meta$seasons$seasonid == years, "description"]
      message(sprintf(
        "No valid years specified, defaulting to this flu season => ID: %s [%s]",
        years, curr_season_descr
      ))
    }
  }

  params$SeasonsDT <- lapply(years, function(i) list(ID = i, Name = as.character(i)))

  tf <- tempfile(fileext = ".zip")
  td <- tempdir()

  on.exit(unlink(tf), TRUE)

  httr::POST(
    url = "https://gis.cdc.gov/grasp/flu2/PostPhase02DataDownload",
    httr::user_agent(.cdcfluview_ua),
    httr::add_headers(
      Origin = "https://gis.cdc.gov",
      Accept = "application/json, text/plain, */*",
      Referer = "https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html"
    ),
    encode = "json",
    body = params,
    # httr::verbose(),
    httr::write_disk(tf)
  ) -> res

  httr::stop_for_status(res)

  nm <- unzip(tf, overwrite = TRUE, exdir = td)

  xdf <- read.csv(nm, skip = 1, stringsAsFactors = FALSE)
  xdf <- .mcga(xdf)

  xdf$weighted_ili <- to_num(xdf$weighted_ili)
  xdf$unweighted_ili <- to_num(xdf$unweighted_ili)
  xdf$age_0_4 <- to_num(xdf$age_0_4)
  xdf$age_25_49 <- to_num(xdf$age_25_49)
  xdf$age_25_64 <- to_num(xdf$age_25_64)
  xdf$age_5_24 <- to_num(xdf$age_5_24)
  xdf$age_50_64 <- to_num(xdf$age_50_64)
  xdf$age_65 <- to_num(xdf$age_65)
  xdf$ilitotal <- to_num(xdf$ilitotal)
  xdf$num_of_providers <- to_num(xdf$num_of_providers)
  xdf$total_patients <- to_num(xdf$total_patients)
  xdf$week_start <- MMWRweek::MMWRweek2Date(xdf$year, xdf$week)

  if (region == "national") xdf$region <- "National"
  if (region == "hhs") xdf$region <- factor(xdf$region, levels = sprintf("Region %s", 1:10))

  class(xdf) <- c("tbl_df", "tbl", "data.frame")

  arrange(suppressMessages(readr::type_convert(xdf)), week_start)

}

utils::globalVariables(c(".", "mmwrid", "season", "seasonid", "week_start", "wk_start", "wk_end", "year_wk_num"))

# CDC U.S. region names to ID map
.region_map <- c(national=3, hhs=1, census=2, state=5)

# CDC hospital surveillance surveillance area name to internal pkg use map
.surv_map <- c(`FluSurv-NET`="flusurv", `EIP`="eip", `IHSP`="ihsp")
.surv_rev_map <- c(flusurv="FluSurv-NET", eip="EIP", ihsp="IHSP")

# CDC P&I mortality GepID mapping
.geoid_map <- c(national="1", state="2", region="3")

# Our bot's user-agent string
.cdcfluview_ua <- "Mozilla/5.0 (compatible; R-cdcvluview Bot/2.0; https://github.com/hrbrmstr/cdcfluview)"

# CDC Basemaps
.national_outline <- "https://gis.cdc.gov/grasp/fluview/FluView2References/Data/US_84.json"
.hhs_subregions_basemap <- "https://gis.cdc.gov/grasp/fluview/FluView2References/Data/HHSRegions_w_SubGroups.json"
.census_divisions_basemap <- "https://gis.cdc.gov/grasp/fluview/FluView2References/Data/CensusDivs_w_SubGroups.json"
.states_basemap <- "https://gis.cdc.gov/grasp/fluview/FluView2References/Data/StatesFluView.json"
.spread_basemap <- "https://gis.cdc.gov/grasp/fluview/FluView8References/Data/States_Territories_labels.json"
.surv_basemap <- "https://gis.cdc.gov/grasp/fluview/FluView1References/data/US_States_w_PR_labels.json"

# CDC Age Groups
.age_grp <- c("0-4 yr", "5-24 yr", "25-64 yr", "65+ yr")

# CDC Virus Groups
.vir_grp <- c("A (Subtyping not Performed)", "A (H1N1)pdm09", "A (Unable to Subtype)",
              "B (Lineage Unspecified)", "A (H1)", "A (H3)", "B (Victoria Lineage)",
              "B (Yamagata Lineage)", "H3N2v")

# Global HTTR timeout
.httr_timeout <- 120

.get_meta <- function() {


  list(
    appversion = jsonlite::unbox("Public"),
    key = jsonlite::unbox(""),
    injson = I(list())
  ) -> body

  httr::POST(
    httr::user_agent(.cdcfluview_ua),
    url = "https://gis.cdc.gov/GRASP/Flu3/PostPhase03DataTool",
    body = jsonlite::toJSON(body),
    encode = "raw",
    httr::accept_json(),
    httr::add_headers(
      `content-type` = "application/json;charset=UTF-8",
      origin = "https://gis.cdc.gov",
      referer = "ttps://gis.cdc.gov/GRASP/Fluview/FluHospRates.html"
    ),
    httr::timeout(.httr_timeout)
  ) -> res

  httr::stop_for_status(res)

  jsonlite::fromJSON(httr::content(res, as = "text"))

}

.mcga <- function(tbl) {

  x <- colnames(tbl)
  x <- tolower(x)
  x <- gsub("[[:punct:][:space:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("(^_|_$)", "", x)
  x <- gsub("^x_", "", x)
  x <- make.unique(x, sep = "_")

  colnames(tbl) <- x

  tbl

}

to_num <- function(x) {
  x <- gsub("%", "", x, fixed=TRUE)
  x <- gsub(">", "", x, fixed=TRUE)
  x <- gsub("<", "", x, fixed=TRUE)
  x <- gsub(",", "", x, fixed=TRUE)
  x <- gsub(" ", "", x, fixed=TRUE)
  suppressWarnings(as.numeric(x))
}

ilinet('hhs', years=2010:2020) %>% write_csv("data/ILInet.csv")
