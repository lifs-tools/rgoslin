test_that("lipid name validation works", {
  expect_equal(rgoslin::isValidLipidName("PC(32:0)"), TRUE)
})

test_that("lipid name parsing works", {
  originalName <- "LPC(34:2;1)"
  df <- rgoslin::parseLipidNames(originalName)
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(df[["Original.Name"]], originalName)
  expect_equal(df[["Normalized.Name"]], "LPC 34:2;O")
  expect_equal(df[["Lipid.Maps.Category"]], "GP")
})

test_that("lipid name parsing with grammar works", {
  originalName <- "LPC(34:1)"
  df <- rgoslin::parseLipidNames(originalName, "SwissLipids")
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(df[["Original.Name"]], originalName)
  expect_equal(df[["Normalized.Name"]], "LPC 34:1")
  expect_equal(df[["Species.Name"]], "LPC 34:1")
  expect_equal(df[["Molecular.Species.Name"]], "NA")
  expect_equal(df[["Sn.Position.Name"]], "NA")
  expect_equal(df[["Structure.Defined.Name"]], "NA")
  expect_equal(df[["Full.Structure.Name"]], "NA")
  expect_equal(df[["Lipid.Maps.Category"]], "GP")
  
  originalName <- "TG(16:1(5E)/18:0/20:2(3Z,6Z))"
  df <- rgoslin::parseLipidNames(originalName, "LipidMaps")
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(df[["Original.Name"]], originalName)
  expect_equal(df[["Normalized.Name"]], "TG 16:1(5E)/18:0/20:2(3Z,6Z)")
  expect_equal(df[["Species.Name"]], "TG 54:3")
  expect_equal(df[["Molecular.Species.Name"]], "TG 16:1_18:0_20:2")
  expect_equal(df[["Sn.Position.Name"]], "TG 16:1/18:0/20:2")
  expect_equal(df[["Structure.Defined.Name"]], "TG 16:1(5)/18:0/20:2(3,6)")
  expect_equal(df[["Full.Structure.Name"]], "TG 16:1(5E)/18:0/20:2(3Z,6Z)")
  expect_equal(df[["Lipid.Maps.Category"]], "GL")
  expect_equal(df[["FA1.DB.Positions"]], "[5E]")
  expect_equal(df[["FA2.DB.Positions"]], "[]")
  expect_equal(df[["FA3.DB.Positions"]], "[3Z, 6Z]")
  
  originalName <- "GalNAcβ1-4(Galβ1-4GlcNAcβ1-3)Galβ1-4Glcβ-Cer(d18:1/24:1(15Z))"
  df <- rgoslin::parseLipidNames(originalName, "LipidMaps")
  # Lipid Maps Species name is currently "Hex(3)-HexNAc(2)-Cer 42:2;O2"
  expect_equal(df[["Species.Name"]], "GalGalGalNAcGlcGlcNAcCer 42:2;O2") 
  expect_equal(df[["Mass"]], 1539.9388, 4)
})

test_that("multiple lipid names parsing works", {
  originalNames <- c("PC(32:0)", "LPC(34:2;1)", "TG(18:1_18:0_16:1)", "TAG 16:1/18:0/20:2")
  df <- rgoslin::parseLipidNames(originalNames)
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(nrow(df), 4)
  expect_equal(as.character(df[1, "Original.Name"]), originalNames[[1]])
  expect_equal(as.character(df[2, "Original.Name"]), originalNames[[2]])
  expect_equal(as.character(df[3, "Original.Name"]), originalNames[[3]])
  expect_equal(as.character(df[4, "Original.Name"]), originalNames[[4]])
  expect_equal(as.character(df[1, "Normalized.Name"]), "PC 32:0")
  expect_equal(as.character(df[2, "Normalized.Name"]), "LPC 34:2;O")
  expect_equal(as.character(df[3, "Normalized.Name"]), "TG 18:1_18:0_16:1")
  expect_equal(as.character(df[4, "Normalized.Name"]), "TG 16:1/18:0/20:2")
})

test_that("multiple lipid names parsing with grammar works", {
  originalNames <- c("PC 32:1","LPC 34:1","TAG 18:1_18:0_16:1")
  df <- rgoslin::parseLipidNames(originalNames, "Goslin")
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(tibble::is_tibble(df), FALSE)
  expect_equal(nrow(df), 3)
  expect_equal(as.character(df[1, "Original.Name"]), originalNames[[1]])
  expect_equal(as.character(df[2, "Original.Name"]), originalNames[[2]])
  expect_equal(as.character(df[3, "Original.Name"]), originalNames[[3]])
  expect_equal(as.character(df[1, "Normalized.Name"]), "PC 32:1")
  expect_equal(as.character(df[2, "Normalized.Name"]), "LPC 34:1")
  expect_equal(as.character(df[3, "Normalized.Name"]), "TG 18:1_18:0_16:1")
})

test_that("lipid name with adduct parsing with grammar works", {
  originalName <- "PC 34:1 [M+H]1+"
  df <- rgoslin::parseLipidNames(originalName, "Goslin")
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(tibble::is_tibble(df), FALSE)
  expect_equal(df[["Original.Name"]], originalName)
  expect_equal(df[["Grammar"]], "Goslin")
  expect_equal(df[["Normalized.Name"]], "PC 34:1")
  expect_equal(df[["Adduct"]], "[M+H]1+")
  expect_equal(df[["Adduct.Charge"]], 1)
  expect_equal(df[["Mass"]], 760.5851, 4)
  expect_equal(df[["Species.Name"]], "PC 34:1")
  expect_equal(df[["Molecular.Species.Name"]], "NA")
  expect_equal(df[["Sn.Position.Name"]], "NA")
  expect_equal(df[["Structure.Defined.Name"]], "NA")
  expect_equal(df[["Full.Structure.Name"]], "NA")
  expect_equal(df[["Lipid.Maps.Category"]], "GP")
  originalName <- "PC 32:1[M+H]+"
  df <- rgoslin::parseLipidNames(originalName, "Goslin")
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(tibble::is_tibble(df), FALSE)
  expect_equal(df[["Original.Name"]], originalName)
  expect_equal(df[["Normalized.Name"]], "PC 32:1")
  expect_equal(df[["Adduct"]], "[M+H]1+")
  expect_equal(df[["Adduct.Charge"]], 1)
  expect_equal(df[["Mass"]], 760.5851, 4)
  expect_equal(df[["Species.Name"]], "PC 32:1")
  expect_equal(df[["Molecular.Species.Name"]], "NA")
  expect_equal(df[["Sn.Position.Name"]], "NA")
  expect_equal(df[["Structure.Defined.Name"]], "NA")
  expect_equal(df[["Full.Structure.Name"]], "NA")
  expect_equal(df[["Lipid.Maps.Category"]], "GP")
})

test_that("parsing many lipid names works", {
  lipidNames <- c(
    "12-HETE",
    "FA(16:0)",
    "BMP 18:1-18:1",
    "LBPA(18:1(11Z)/0:0/32:5(14Z,17Z,20Z,23Z,26Z)/0:0)",
    "CDPDAG 18:1-18:1",
    "Cer 18:1;2/16:0",
    "Cer(d18:1/18:0)",
    "CL(16:0/16:0/16:0/20:0)",
    "CL 18:3(9Z,12Z,15Z)/16:0/22:5(16Z,10Z,19Z,13Z,7Z)/20:0",
    "DG(18:2_20:4)",
    "DGDG 16:0-16:1",
    "GB3 18:1;2/24:1",
    "Gb3(d18:1(4E)/24:1(15Z))",
    "Hex2Cer 18:1;2/12:0",
    "HexCer(d18:1/20:0)",
    "GalCer(d18:1(4E)/20:0)",
    "GlcCer(d18:1(4E)/20:0)",
    "LCB 17:1;2",
    "LPC(20:3)",
    "LPC(O-22:1)",
    "LPE O 19:1p",
    "MAG 16:0",
    "MLCL 18:1-18:1-18:1",
    "MLCL (20:2(11Z,14Z)/20:1(11Z)/20:0/0:0)",
    "M(IP)2C(t20:0/26:0(2OH))",
    "Palmitic acid",
    "PC(18:2_20:4)",
    "PC(O-40:7)",
    "PC(P-30:0)",
    "PE 16:2-18:3;1",
    "PE 16:2/18:3;1",
    "PE 18:3;1-16:2",
    "PE O 18:0a/22:6",
    "PE(P-16:0/22:6)",
    "PEt 16:0-18:1",
    "PS 18:2-22:1",
    "PIP2 21:0-22:6",
    "SHexCer 18:0;3/26:0;1",
    "SM(d35:1)",
    "SM d35:1 [M+H]1+",
    "TG(14:0_16:0_18:1)",
    "PC(21:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z))",
    "PE(16:2(9Z,12Z)/18:1(6Z))",
    "PIP[4'] 6:0/22:5(16Z,10Z,19Z,13Z,7Z)",
    "TG(18:1/10:0/16:0)",
    "TG(20:0/22:3(10Z,13Z,16Z)/22:5(7Z,10Z,13Z,16Z,19Z))[iso6]"
  )
  for(i in 1:length(lipidNames)) {
    df <- rgoslin::parseLipidNames(lipidNames[[i]])[1,]
    expect_equal(is.data.frame(df), TRUE)
    expect_equal(tibble::is_tibble(df), FALSE)
    expect_equal(as.character(df[1, "Original.Name"]), lipidNames[[i]])
    expect_true(!is.na(df[1, "Normalized.Name"]))
  }
  df2 <- rgoslin::parseLipidNames(lipidNames)
  expect_equal(46, nrow(df2))
  expect_equal(is.data.frame(df2), TRUE)
  expect_equal(tibble::is_tibble(df2), FALSE)
  expect_equal(0, sum(is.na(df2[, "Message"])))
  expect_equal(0, sum(is.na(df2[, "Normalized.Name"])))
})

test_that("getting the list of supported parsers/grammars works", {
  expectedGrammars <- c(
    "FattyAcids",
    "Goslin",
    "HMDB",
    "LipidMaps",
    "Shorthand2020",
    "SwissLipids"
  )
  availableGrammars <- rgoslin::listAvailableGrammars() %>% sort()
  expect_setequal(expectedGrammars, availableGrammars)
})

test_that("lipid level works", {
  l <- rgoslin::parseLipidNames("PE 16:1(6Z)/16:0;5OH[R],8OH;3oxo", "Shorthand2020")[1,]
  expect_equal("COMPLETE_STRUCTURE", l[["Level"]])
  
  l <- rgoslin::parseLipidNames("PE 16:1(6Z)/16:0;5OH,8OH;3oxo", "Shorthand2020")[1,]
  expect_equal("FULL_STRUCTURE", l[["Level"]])
  
  l <- rgoslin::parseLipidNames("PE 16:1(6)/16:0;(OH)2;oxo", "Shorthand2020")[1,]
  expect_equal("STRUCTURE_DEFINED", l[["Level"]])
  
  l <- rgoslin::parseLipidNames("PE 16:1/16:1;O3", "Shorthand2020")[1,]
  expect_equal("SN_POSITION", l[["Level"]])
  
  l <- rgoslin::parseLipidNames("PE 16:1_16:1;O3", "Shorthand2020")[1,]
  expect_equal("MOLECULAR_SPECIES", l[["Level"]])
  
  l <- rgoslin::parseLipidNames("PE 32:1;O3", "Shorthand2020")[1,]
  expect_equal("SPECIES", l[["Level"]])
})

test_that("cyclopropane works", {
  l <- rgoslin::parseLipidNames("FA 19:0;[11-13cy3:0]", "Shorthand2020")[1,]
  expect_equal("FULL_STRUCTURE", l[["Level"]])
})

test_that("DB count for Fa is correct", {
  l <- rgoslin::parseLipidNames("CAR 18:1")[1,]
  expect_equal(1, l[["FA1.DB"]])
  expect_equal(425.3505, l[["Mass"]], tolerance <- 1e-06)
})

test_that("LCB and FAs are distinguished", {
  l <- rgoslin::parseLipidNames("Cer d18:1/24:0")[1,]
  expect_equal(1, l[["LCB.Position"]])
  expect_equal(18, l[["LCB.C"]])
  expect_equal(1, l[["LCB.DB"]])
  expect_equal(2, l[["LCB.OH"]])
  expect_equal("LCB", l[["LCB.Bond.Type"]])
  expect_equal(2, l[["FA1.Position"]])
  expect_equal(24, l[["FA1.C"]])
  expect_equal(0, l[["FA1.DB"]])
})

test_that("Hydroxyl group counts are proper", {
  l <- rgoslin::parseLipidNames("Cer 36:1;2")[1,]
  expect_equal(2, l[["Total.OH"]])
  l <- rgoslin::parseLipidNames("Cer d36:1")[1,]
  expect_equal(2, l[["Total.OH"]])
  l <- rgoslin::parseLipidNames("Cer 18:1;2/18:0")[1,]
  expect_equal(2, l[["Total.OH"]])
  l <- rgoslin::parseLipidNames("Cer d18:1/18:0")[1,]
  expect_equal(2, l[["Total.OH"]])
  l <- rgoslin::parseLipidNames("Cer 18:1;(OH)2/18:0")[1,]
  expect_equal(2, l[["Total.OH"]])
})

test_that("IUPAC Fatty Acid Names are parsed correct", {
  l <- rgoslin::parseLipidNames("5-methyl-octadecanoic acid")[1,]
  expect_equal("FA 18:0;5Me", l[["Normalized.Name"]])
  l <- rgoslin::parseLipidNames("2-docosyl-3-hydroxy-28,29-epoxy-30-methyl-pentacontanoic acid")[1,]
  expect_equal("FA 50:0;2(22:0);28Ep;30Me;3OH", l[["Normalized.Name"]])
  l <- rgoslin::parseLipidNames("11R-hydroxy-9,15-dioxo-2,3,4,5-tetranor-prostan-1,20-dioic acid")[1,]
  expect_equal("FA 15:0;15COOH;[4-8cy5:0;7OH;5oxo];11oxo", l[["Normalized.Name"]])
  l <- rgoslin::parseLipidNames("N-((+/-)-8,9-dihydroxy-5Z,11Z,14Z-eicosatrienoyl)-ethanolamine", "FattyAcids")[1,]
  expect_equal("NAE 20:3(5Z,11Z,14Z);8OH,9OH", l[["Normalized.Name"]])
  l <- rgoslin::parseLipidNames("N-ethyl-5Z,8Z,11Z,14Z-eicosatetraenoyl amine", "FattyAcids")[1,]
  expect_equal("NA 2:0/20:4(5Z,8Z,11Z,14Z)", l[["Normalized.Name"]])
  expect_equal("NA 2:0/20:4(5,8,11,14)", l[["Structure.Defined.Name"]])
  expect_equal("NA 2:0/20:4", l[["Sn.Position.Name"]])
  l <- rgoslin::parseLipidNames(l[["Sn.Position.Name"]], "Shorthand2020")[1,]
  expect_equal("NA 2:0_20:4", l[["Molecular.Species.Name"]])
  expect_equal("NA 22:4", l[["Species.Name"]])
})

test_that("Sterols work", {
  l <- rgoslin::parseLipidNames("Desmosterol", "Goslin")
  expect_equal("ST 27:2;O", l[["Normalized.Name"]])
})

test_that("isValidLipidName creates warning on invalid name",{
  expect_warning(rgoslin::isValidLipidName("PX 40:1"), "Parsing of lipid name 'PX 40:1' caused an exception: Lipid not found")
})

test_that("isValidLipidName stops on invalid input",{
  expect_error(rgoslin::isValidLipidName(c(1)),"'lipidName' must be a string")
})

test_that("weird input creates messages", {
  expect_message(rgoslin::parseLipidNames(c("A")), "Encountered an error while parsing 'A': Expecting a single string value: ")
  expect_message(rgoslin::parseLipidNames(c("")), "Encountered an error while parsing '': Expecting a single string value: ")
  expect_message(rgoslin::parseLipidNames(c(""), "LipidMaps"), "Encountered an error while parsing '' with grammar 'LipidMaps': Expecting a single string value: ")
  # try to parse a) with wrong grammar, b) with invalid input
  l <- rgoslin::parseLipidNames(c("Cer 36:1;2", 5), "Shorthand2020")
  expect_equal(nrow(l), 2)
  expect_true(is.na(l[1,"Normalized.Name"]))
  expect_true(is.na(l[2,"Normalized.Name"]))
  expect_equal("NOT_PARSEABLE", l[1,"Grammar"])
  expect_equal("NOT_PARSEABLE", l[2,"Grammar"])
  l2 <- rgoslin::parseLipidNames(c("Cer 36:1;2", 5))
  expect_equal(nrow(l2), 2)
  expect_equal("NOT_PARSEABLE", l[2,"Grammar"])
  expect_true(is.na(l[2,"Normalized.Name"]))
  testthat::expect_error(rgoslin::parseLipidNames(c(1,2)), "lipidNames must not contain numbers only!")
})
