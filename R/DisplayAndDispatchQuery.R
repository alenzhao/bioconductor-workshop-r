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

#' Display and dispatch a query to BigQuery.
#'
#' This convenience method will read SQL from a local file or a url.  It will
#' then perform the text substitutions, if any, and dispatch the query to the cloud.
#'
#' @param queryUri File path or url to file containing SQL code.
#' @param project The ID of a Google Cloud Platform project of which the user is a member.
#' @param replacements A list of key/value pairs where the elements are the "values", and
#'  the name of the element is the corresponding key i.e. replacements[[KEY]] == VALUE.
#'  For each pair, if the KEY is found in the SQL using an exact string match criteria,
#'  it will be replaced with the VALUE.
#' @return The dataframe of query results.
#' @export
DisplayAndDispatchQuery <- function(queryUri, project, replacements=list()) {
  if (missing(queryUri)) {
    stop("Pass the file path or url to the file containing the query.")
  }
  if(missing(project)) {
    stop("Pass the project id of your Google Cloud Platform project.")
  }

  if (grepl("^https.*", queryUri)) {
    # Read the query from a remote location.
    querySql <- RCurl::getURL(queryUri, ssl.verifypeer=FALSE)
  } else {
    # Read the query from the local filesystem.
    querySql <- readChar(queryUri, nchars=1e6)
  }

  # If applicable, substitute values in the query template.
  for(replacement in names(replacements)) {
    querySql <- gsub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }

  # Display the query to the terminal.
  cat(querySql)

  # Dispatch the query to BigQuery.
  bigrquery::query_exec(querySql, project)
}
