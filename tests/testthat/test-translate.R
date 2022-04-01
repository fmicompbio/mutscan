test_that("translateString works", {
  code <- Biostrings::GENETIC_CODE
  for (i in seq_along(code)) {
    expect_identical(translateString(names(code)[i]), unname(code[i]))
  }
  expect_identical(translateString(paste(names(code), collapse = "_")),
                   paste(code, collapse = "_"))
  expect_identical(translateString(paste(names(code), collapse = "")),
                   paste(code, collapse = ""))
  expect_identical(translateString(tolower(paste(names(code), collapse = ""))),
                   paste(code, collapse = ""))
  expect_identical(translateString(paste(gsub("T", "U", names(code)), collapse = "")),
                   paste(code, collapse = ""))
  expect_identical(translateString("AAN,AAA.A"), "XXX")
})
