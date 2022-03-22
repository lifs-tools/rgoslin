#' @useDynLib rgoslin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Check lipid name.
#'
#' \code{isValidLipidName} checks the provided lipid name against the built-in
#' grammars.
#' Will return FALSE if none of the parsers was able to parse the provided name
#' successfully. Will stop execution via \code{stop} if non character input is 
#' detected.
#' @param lipidName The lipid name to check.
#' @examples
#' isValidLipidName("PC 32:1")
#' isValidLipidName("PC(32:1)")
#' isValidLipidName("PCX(32:1)")
#' @return TRUE if the lipidName could be parsed, FALSE otherwise.
#' @export
isValidLipidName <- function(lipidName) {
    if(!is.character(lipidName)) {
        stop("'lipidName' must be a string")
    }
    tryCatch({
        return(rcpp_is_valid_lipid_name(lipidName))
    })
}

#' Parse multiple lipid names and return a data frame with the results.
#'
#' \code{parseLipidNames} reads the provided lipid names vector and returns
#' structural information as a data frame. Will return a cell with the
#' "Grammar" column set to "NOT_PARSEABLE" if none of the parsers was able to
#' parse the provided name successfully. Will stop execution via \code{stop} if 
#' invalid non character input is detected or fatal errors are encountered 
#' during parsing.
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
    if (is.numeric(lipidNames)) {
        stop("lipidNames must not contain numbers only!")
    }
    for (i in seq_along(lipidNames)) {
        df <- tryCatch({
            if (is.null(grammar)) {
                as.data.frame(
                    rcpp_parse_lipid_name(
                        as.character(lipidNames[[i]])
                    )
                )
            } else {
                as.data.frame(
                    rcpp_parse_lipid_name_with_grammar(
                        as.character(lipidNames[[i]]), grammar
                    )
                )
            }
        }, warning = function(warn) {
            grammarStr <- ifelse(
                is.null(grammar),
                "",
                paste0(
                    " with grammar '",
                    grammar,
                    "'"
                )
            )
            message(
                "Encountered a warning while parsing '",
                lipidNames[[i]],
                "'",
                grammarStr,
                ": ",
                warn$message
            )
            data.frame(
                "Normalized.Name" = NA,
                "Original.Name" = as.character(lipidNames[[i]]),
                "Grammar" = "NOT_PARSEABLE",
                "Message" = warn$message
            )
        }, error = function(err) {
            grammarStr <- ifelse(
                is.null(grammar),
                "",
                paste0(
                    " with grammar '",
                    grammar,
                    "'"
                )
            )
            message(
                "Encountered an error while parsing '",
                lipidNames[[i]],
                "'",
                grammarStr,
                ": ",
                err$message
            )
            data.frame(
                "Normalized.Name" = NA,
                "Original.Name" = as.character(lipidNames[[i]]),
                "Grammar" = "NOT_PARSEABLE",
                "Message" = err$message
            )
        })
        namesList[[i]] <- df
    }
    return(dplyr::bind_rows(namesList))
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
