
# surveybootstrap 0.0.3

* fixed some dependency errors that suddenly arose (explicitly use plyr:: in a few places)
* dropped now-deprecated group_indices_() in favor of cur_group_index()
* changed docs to use Markdown roxygen
* added MU284 dataset for examples and tests
* added examples to exported function docs

# surveybootstrap 0.0.2

* moved unit tests for goc and helper functions in from networkreporting
* added rescaled.bootstrap.weights helper function for case when a tibble with all of the rescaled bootstrap weights is useful
* added dependencies to tidyverse, magrittr, and purrr

# surveybootstrap 0.0.1

* first CRAN version, May 2016
