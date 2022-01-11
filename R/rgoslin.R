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
    }, error = function(err) {
        message(
            "Could not parse ",
            lipidName,
            " with any of the available parsers!"
        )
        return(FALSE)
    })
}

#' Parse lipid name.
#'
#' \code{parseLipidName} reads the provided lipid name and returns structural
#' information as a named vector. Will return a vector with the "Grammar"
#' element set to "NOT_PARSEABLE" if none of the parsers was able to parse the
#' provided name successfully. If any error was raised, returns an empty vector.
#' @param lipidName The lipid name to parse.
#' @examples
#' parseLipidName("PC 32:1")
#' parseLipidName("LPC 34:1")
#' parseLipidName("TG(18:1_18:0_16:1)")
#' @return Data frame with details of the lipid, empty data frame otherwise.
#' @export
parseLipidName <- function(lipidName) {
    tryCatch({
        return(as.data.frame(rcpp_parse_lipid_name(lipidName)))
    }, error = function(err) {
        message(
            "Could not parse ",
            lipidName,
            " with any of the available parsers!"
        )
        return(data.frame())
    })
}

#' Parse lipid name with a specific grammar.
#'
#' \code{parseLipidName} reads the provided lipid name based on grammar and
#' returns structural information as a named vector. Will return a vector with
#' the "Grammar" element set to "NOT_PARSEABLE" if none of the parsers was able
#' to parse the provided name with the provided grammar successfully.If any
#' error was raised, returns an empty vector.
#' @param lipidName The lipid name to parse.
#' @param grammar The grammar to use. One of "Goslin", "GoslinFragments",
#' "SwissLipids", "LipidMaps", "HMDB"
#' @examples
#' parseLipidNameWithGrammar("PC 32:1", "Goslin")
#' parseLipidNameWithGrammar("LPC(34:1)", "SwissLipids")
#' parseLipidNameWithGrammar("TG(18:1_18:0_16:1)", "LipidMaps")
#' @return character vector with details of the lipid, empty vector otherwise.
#' @export
parseLipidNameWithGrammar <- function(lipidName, grammar) {
    tryCatch({
        return(as.data.frame(
            rcpp_parse_lipid_name_with_grammar(lipidName, grammar)
        ))
    }, error = function(err) {
        message(
            "Could not parse ",
            lipidName,
            " with any of the available parsers!"
        )
        return(data.frame())
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
#' @examples
#' parseLipidNames(c("PC 32:1","LPC 34:1","TG(18:1_18:0_16:1)"))
#' @return Data frame where each row reports the parsing result of each element
#' in lipidNames.
#' @export
parseLipidNames <- function(lipidNames) {
    tryCatch({
        namesList <- list()
        for (i in seq_along(lipidNames)) {
            tryCatch({
                namesList[[i]] <-
                    rcpp_parse_lipid_name(as.character(lipidNames[[i]]))
            }, error = function(err) {
                message(
                    "Could not parse ",
                    lipidNames[[i]],
                    " with any of the available parsers!"
                )
            })
        }
        return(dplyr::bind_rows(namesList))
    }, error = function(err) {
        msg <- paste(
            "Could not parse the provided lipid names",
            paste0(lipidNames, collapse=","),
            "with any of the available parsers!",
        )
        message(
            msg
        )
        return(data.frame())
    })
}

#' Parse multiple lipid names with the provided grammar and return a data frame
#' with the results.
#'
#' \code{parseLipidNamesWithGrammar} reads the provided lipid names vector and
#' parses it against the given grammar. It returns structural information as a
#' data frame. Will return a cell with the "Grammar" column set to
#' "NOT_PARSEABLE" if none of the parsers was able to parse the provided name
#' successfully. If any error was raised, returns an empty data frame.
#' @param lipidNames The vector of lipid names to parse.
#' @param grammar The grammar to use. One of "Goslin", "GoslinFragments",
#' "SwissLipids", "LipidMaps", "HMDB".
#' @examples
#' parseLipidNamesWithGrammar(c("PC 32:1","LPC 34:1","TG(18:1_18:0_16:1)"),
#' "Goslin")
#' @return Data frame where each row reports the parsing result of each element
#' in lipidNames.
#' @export
parseLipidNamesWithGrammar <- function(lipidNames, grammar) {
    tryCatch({
        namesList <- list()
        for (i in seq_along(lipidNames)) {
            tryCatch({
                namesList[[i]] <-
                    as.data.frame(
                        rcpp_parse_lipid_name_with_grammar(
                            lipidNames[[i]], grammar
                        )
                    )
            }, error = function(err) {
                message(
                    "Could not parse ",
                    lipidNames[[i]],
                    " with any of the available parsers!"
                )
            })
        }
        return(do.call(rbind, namesList))
    }, error = function(err) {
        msg <- paste(
            "Could not parse the provided lipid names",
            paste0(lipidNames, collapse=","),
            "with any of the available parsers!",
        )
        message(
            msg
        )
        return(data.frame())
    })
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
    tryCatch({
        return(rcpp_list_available_grammars())
    }, error = function(err) {
        message("Could not access supported grammar names!")
        return(c())
    })
}
