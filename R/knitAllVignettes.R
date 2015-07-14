# Copyright 2015 Google Inc. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' These vignettes require an oauth2 flow and therefore cannot be built as
#' part of the normal package build mechanism.
#'
#' Instead we build them manually and check in the rendered results.
knitAllVignettes <- function() {
  require(knitr)
  if(FALSE == grepl('doc$', getwd())) {
    stop("be sure to setwd('PATH/TO/inst/doc') before running this command.")
  }
  lapply(c("BigQueryDemo.Rmd", "AllModalitiesDemo.Rmd"), function(rmd) {
    knit(rmd)
    knit2html(rmd)
    purl(rmd, documentation=2)
  }
  )
}
