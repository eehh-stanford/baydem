## Download and read the IntCal13, Marine13, and SHCal13 calibration curves
library(tidyverse)
library(magrittr)
library(httr)
library(rvest)
library(writexl)

# IntCal13
download.file("http://www.radiocarbon.org/IntCal13%20files/intcal13.14c", "data-raw/intcal13.14c")
intcal13 <- readr::read_csv("data-raw/intcal13.14c", skip = 9)[-1, ] %>%
  dplyr::rename(`CAL BP` = `# CAL BP`) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `CAL BP`,
      `14C age`,
      Error
    ),
    as.integer
  ) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `Delta 14C`,
      Sigma
    ),
    as.numeric
  )

usethis::use_data(intcal13,
  overwrite = T
)

# Marine13
download.file("http://www.radiocarbon.org/IntCal13%20files/marine13.14c", "data-raw/marine13.14c")
marine13 <- readr::read_csv("data-raw/marine13.14c", skip = 9)[-1, ] %>%
  dplyr::rename(`CAL BP` = `# CAL BP`) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `CAL BP`,
      `14C age`,
      Error
    ),
    as.integer
  ) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `Delta 14C`,
      Sigma
    ),
    as.numeric
  )

usethis::use_data(marine13,
  overwrite = T
)

# SHCal13
download.file("http://www.radiocarbon.org/IntCal13%20files/shcal13.14c", "data-raw/shcal13.14c")
shcal13 <- readr::read_csv("data-raw/shcal13.14c", skip = 9)[-1, ] %>%
  dplyr::rename(`CAL BP` = `# CAL BP`) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `CAL BP`,
      `14C age`,
      Error
    ),
    as.integer
  ) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `Delta 14C`,
      Sigma
    ),
    as.numeric
  )

usethis::use_data(shcal13,
  overwrite = T
)

# Raw data for IntCal13
intcal13_raw <- c("cooked", "raw", "sets", "refs") %>%
  magrittr::set_names(., .) %>%
  purrr::map(function(x) {
    httr::POST("http://intcal.qub.ac.uk/intcal13/query/query.php",
      body = list(query = paste0("select * from ", x))
    ) %>%
      httr::content() %>%
      rvest::html_nodes("table") %>%
      .[2] %>%
      rvest::html_table(fill = TRUE, header = TRUE) %>%
      magrittr::extract2(1) %>%
      tibble::as_tibble()
  }) %T>%
  writexl::write_xlsx("data-raw/intcal13_raw.xlsx")

usethis::use_data(intcal13_raw,
  overwrite = T
)

# Raw data for Marine13
marine13_raw <- c("details", "timedep", "taxa", "class", "refs") %>%
  magrittr::set_names(., .) %>%
  purrr::map(function(x) {
    httr::POST("http://intcal.qub.ac.uk/marine/query/query.php",
      body = list(query = paste0("select * from ", x))
    ) %>%
      httr::content() %>%
      rvest::html_nodes("table") %>%
      .[2] %>%
      rvest::html_table(fill = TRUE, header = TRUE) %>%
      magrittr::extract2(1) %>%
      tibble::as_tibble()
  }) %T>%
  writexl::write_xlsx("data-raw/marine13_raw.xlsx")

usethis::use_data(marine13_raw,
  overwrite = T
)

# Raw data for SHCal13
shcal13_raw <- c("cooked", "raw", "sets") %>%
  magrittr::set_names(., .) %>%
  purrr::map(function(x) {
    httr::POST("http://intcal.qub.ac.uk/shcal13/query/query.php",
      body = list(query = paste0("select * from ", x))
    ) %>%
      httr::content() %>%
      rvest::html_nodes("table") %>%
      .[2] %>%
      rvest::html_table(fill = TRUE, header = TRUE) %>%
      magrittr::extract2(1) %>%
      tibble::as_tibble()
  }) %T>%
  writexl::write_xlsx("data-raw/shcal13_raw.xlsx")

usethis::use_data(shcal13_raw,
  overwrite = T
)
