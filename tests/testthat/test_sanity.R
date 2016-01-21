############
# some basic sanity checks
############
library(miRNAtap)


context('main function high level tests')


test_that("number of output cols corresponds to input parameters", {
    expect_equal(dim(getPredictedTargets('let-7a',species='mmu',
                                        method='geom'))[2],6)
    expect_equal(dim(getPredictedTargets('let-7a',species='mmu', method='geom',
                                        sources=c('pictar','diana',
                                                    'targetscan')
                                        ))[2], 5)
    expect_equal(dim(getPredictedTargets('let-7a',species='mmu', method='geom',
                                        sources=c('pictar','diana'),
                                        min_src = 1))[2], 4)
    expect_equal(dim(getPredictedTargets('let-7a',species='mmu', method='geom',
                                        sources=c('pictar'),
                                        min_src = 1))[2], 3)
})


test_that("stupid parameters return null", {
    expect_null(getPredictedTargets('notamirn-4',species='mmu', method='geom'))
    expect_null(getPredictedTargets('let-7a',species='mmu', method='geom',
                                    sources=c('pictar'), min_src = 2))
    expect_null(getPredictedTargets('let-7a',species='mmu', method='geom',
                                    sources=c('not_a_source'), min_src = 1))
})

test_that("increasing min_src decreases recall", {
    res1 <- getPredictedTargets("let-7a", species = "mmu",
                                method = "geom", min_src = 1)
    res2 <- getPredictedTargets("let-7a", species = "mmu",
                                method = "geom", min_src = 2)
    res3 <- getPredictedTargets("let-7a", species = "mmu",
                                method = "geom", min_src = 3)
    res4 <- getPredictedTargets("let-7a", species = "mmu",
                                method = "geom", min_src = 4)
    
    expect_more_than(dim(res1)[1],
                    dim(res2)[1]-1)
    expect_more_than(dim(res2)[1],
                    dim(res3)[1]-1)
    expect_more_than(dim(res3)[1],
                    dim(res4)[1]-1)
    expect_more_than(dim(res4)[1],
                    0)
})