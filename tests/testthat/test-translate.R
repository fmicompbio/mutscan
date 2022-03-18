test_that("translate works", {
  code <- Biostrings::GENETIC_CODE
  for (i in seq_along(code)) {
    expect_identical(translate(names(code)[i]), unname(code[i]))
  }
  expect_identical(translate(paste(names(code), collapse = "_")),
                   paste(code, collapse = "_"))
  expect_identical(translate(paste(names(code), collapse = "")),
                   paste(code, collapse = ""))
  expect_identical(translate(tolower(paste(names(code), collapse = ""))),
                   paste(code, collapse = ""))
  expect_identical(translate(paste(gsub("T", "U", names(code)), collapse = "")),
                   paste(code, collapse = ""))
  expect_identical(translate("AAN,AAA.A"), "XXX")
})
