#' @useDynLib rgoslin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Check lipid name.
#'
#' \code{isValidLipidName} checks the provided lipid name against the built-in
#' grammars.
#' Will return FALSE if none of the parsers was able to parse the provided name
#' successfully or if any error was raised.
#' @param lipidName The lipid name to check.
#' @examples
#' isValidLipidName("PC 32:1")
#' isValidLipidName("PC(32:1)")
#' isValidLipidName("PCX(32:1)")
#' @return TRUE if the lipidName could be parsed, FALSE otherwise.
#' @export
isValidLipidName <- function(lipidName) {
    tryCatch({
        return(rcpp_is_valid_lipid_name(lipidName))
    })
}

#' Parse multiple lipid names and return a data frame with the results.
#'
#' \code{parseLipidNames} reads the provided lipid names vector and returns
#' structural information as a data frame. Will return a cell with the
#' "Grammar" column set to "NOT_PARSEABLE" if none of the parsers was able to
#' parse the provided name successfully. If any error was raised, returns an
#' empty data frame.
#' @param lipidNames The vector of lipid names to parse.
#' @param grammar The grammar to use. One of "Goslin", "GoslinFragments",
#' "SwissLipids", "LipidMaps", "HMDB", "FattyAcids". Call
#' \code{listAvailableGrammars()} for a complete list of available grammars.
#' If \code{grammar} is omitted or \code{NULL} is passed as a parameter, all
#' available grammars / parsers will be tested. The first successful one will
#' win. If all parsers fail, the "Messages" column in the returned data frame
#' will contain the last parsers message.
#' @examples
#' parseLipidNames(c("PC 32:1","LPC 34:1","TG(18:1_18:0_16:1)"))
#' parseLipidNames(c("Cer(d18:1(8Z)/24:0)", grammar = "LipidMaps"))
#' @return Data frame where each row reports the parsing result of each element
#' in lipidNames.
#' @export
parseLipidNames <- function(lipidNames, grammar = NULL) {
    namesList <- list()
    for (i in seq_along(lipidNames)) {
        if (is.null(grammar)) {
            tryCatch({
            namesList[[i]] <-
                as.data.frame(
                    rcpp_parse_lipid_name(
                        as.character(lipidNames[[i]])
                    )
                )
            })
        } else {
            tryCatch({
                namesList[[i]] <-
                    as.data.frame(
                        rcpp_parse_lipid_name_with_grammar(
                            as.character(lipidNames[[i]]), grammar
                        )
                    )
            })
        }
    }
    return(do.call(rbind, namesList))
}

#' Return the list of grammars supported by goslin.
#'
#' \code{listAvailableGrammars} returns the list of grammars that the
#' underlying cppgoslin library supports.
#' @examples
#' listAvailableGrammars()
#' @return the list of grammars
#' @export
listAvailableGrammars <- function() {
    return(rcpp_list_available_grammars())
}
