# rgosing 1.5.0

Please note that this Bioconductor version is based on Goslin version 2.0.0.
See the [Goslin repository](https://github.com/lifs-tools/goslin) and [Goslin C++ repository](https://github.com/lifs-tools/cppgoslin) for more details.

## BioConductor 3.18 - Changes in 1.5.0

### Improvements

- Better handling of mediators
- Translating gangliosides into new nomenclature structure
- Updated HMDB grammar for parsing mediators in lipid names: e.g., PA(P-16:0/LTE4)
- Parsing of adducts with heavy labeled isotopes now possible

### Bug Fixes

- Minor bug fixes

# rgoslin 1.4.0 

Please note that this Bioconductor version is based on Goslin version 2.0.0.
See the [Goslin repository](https://github.com/lifs-tools/goslin) and [Goslin C++ repository](https://github.com/lifs-tools/cppgoslin) for more details.

## BioConductor 3.17 - Changes in 1.4.0

### Improvements

- Improved handling of Glycosphingolipids and carbohydrates
- Improved headgroup normalization for Glycosphingolipids.
- Added PMeOH.
- Added TG-EST (estolide) Estolides [GL0305].
- Added more sterol variants.

### Bug Fixes

- Fixed classification of SB1a as Sulfoglycosphingolipids (sulfatides) [SP0602].
- Fixed classification of SHex2Cer as Sulfoglycosphingolipids (sulfatides) [SP0602].
- Fixed classification of SMGDG as Glycosylalkylacylglycerols [GL0502], added synonym seminolipid.
- Fixed classification of SQDG Glycosyldiacylglycerols [GL0501].
- Fixed classification of sterols.

## BioConductor 3.16 - Changes in 1.2.0

- No noteworthy changes.

## BioConductor 3.15 - Changes in 1.0.0

### Improvements

- Reduced memory consumption.
- Added 'ChE' abbreviation.
- Added FG hydroperoxy to mediator nomenclature, refinement of mediators.
- Added more sphingosine and sphinganine synonyms.
- Added more ether dialects to LipidMaps grammar.
- Improved handling for SP without explicit OH description.
- Added Sa So support.
- Updated old SP shortcuts.
- Added CholE as abbreviation for cholesterol esters.
- Modifications and improvements for Windows.
- Added column of elements to functional group list and class.
- Added 'ChoE'.
- Added functional group butylperoxy -> BOO.

### Bug Fixes

- Fixed handling of LIPID MAPS SP notation.
- Fixed critical bug when parsing LIPID MAPS names.
- Fixed implicit hydroxy count.
- Fixed ACer rule for species level.
- Fixed lcb rule in LipidMaps grammar.
- Fixed S1P and Sa1P handling.
- Fixed gangliosides in Goslin grammar.
- Fixed correct handling of dummy FAs during sorting.
- Fixed segmentation fault in FA parser event handler.

## Changes in 0.99.1

- The column names within the data frames returned from the `parse*` methods now use column names with dots instead of spaces. This makes it easier to use the column names unquoted within other R expressions.
- All `parse*` methods now return data frames.
- The `Messages` column has been added to capture parser messages. If parsing succeeds, this will contain `NA` and `Normalized.Name` will contain the normalized lipid shorthand name.
- Parser implementations have been updated to reflect the latest lipid shorthand nomenclature changes. Please see the [Goslin repository](https://github.com/lifs-tools/goslin) for more details.
- Exceptions in the C++ part of the library are captured as warnings in R. However, if you parse multiple lipid names, exceptions will not stop the parsing process.

