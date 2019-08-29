ternSimNoHerr <- predOutn %>%
  + filter(Yr>99)
ternSim <-ggplot(ternSimNoHerr, aes(x=Yr, y=PredN)) + geom_point()
ternSim
mean(ternSimNoHerr$PredN)
lm(PredN ~ Yr, ternSimNoHerr)

#real tern pop
mean(ternpopproddiet$productivity[ternpopproddiet$species=="Common"], na.rm=T)
#[1] 1.018827