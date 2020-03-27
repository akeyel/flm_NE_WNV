test_that('Code will run with package sample data', {
  
  human.data = weather.data = NA
  start.year = 2002
  flm.inputs = list(flm::NE_county_pops, flm::spi, flm::spei, flm::sampledat, flm::NEdat, start.year)
  weekinquestion = "2018-02-01" #**# Check date formatting, and what date format comes in from DFMIP
  results.path = '/testing_temporary'
  analysis.districts = unique(flm::sampledat$County)
  n.draws = 1
  
  flm.out = run.flm(human.data, weather.data, flm.inputs, weekinquestion, results.path, in.seed = 4872957)
  flm.results = flm.out[[1]]
  expect_equal(length(flm.out), 3)
  expect_equal(nrow(flm.results), 93)
  expect_equal(is.na(flm.results$predcases[1]), TRUE) #**# For now, this is giving me an error. This test is here to remind me to change the unit test once the NaN is fixed.
  
  #**# Can add additional unit tests as desired
  
  # Remove temporary testing path
  unlink(results.path)
})


