# rgoslin 0.99

Please note that this Bioconductor version is based on Goslin version 2.0.0.
See the [Goslin repository](https://github.com/lifs-tools/goslin) for more details.

## Changes in 0.99.1

- The column names within the data frames returned from the `parse*` methods now use column names with dots instead of spaces. This makes it easier to use the column names unquoted within other R expressions.
- All `parse*` methods now return data frames.
- The `Messages` column has been added to capture parser messages. If parsing succeeds, this will contain `NA` and `Normalized.Name` will contain the normalized lipid shorthand name.
- Parser implementations have been updated to reflect the latest lipid shorthand nomenclature changes. Please see the [Goslin repository](https://github.com/lifs-tools/goslin) for more details.
- Exceptions in the C++ part of the library are captured as warnings in R. However, if you parse multiple lipid names, exceptions will not stop the parsing process.

