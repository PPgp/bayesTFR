e <- new.env()
load('iso3166.rda', envir = e)
iso3166 <- e$iso3166
iso3166 <- iso3166[,c(1,3,4,5)]
colnames(iso3166) <- c('name', 'charcode', 'charcode3', 'uncode')
rm(e)
