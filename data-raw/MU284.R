# NOTE: the file mu284-surveys.r in this directory has the code
# that was originally used to generate these data

load('data-raw/MU284.RData')

usethis::use_data(MU284, overwrite = TRUE)
usethis::use_data(MU284.surveys, overwrite = TRUE)
