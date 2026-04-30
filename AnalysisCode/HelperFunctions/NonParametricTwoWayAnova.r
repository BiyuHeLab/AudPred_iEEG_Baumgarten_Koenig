library(ARTool)

df <- read.csv("FTPLrating_for_ART.csv")

df$Subj <- factor(df$Subj)
df$Predp34 <- factor(df$Predp34)
df$p34 <- factor(df$p34)
df$ToneDur <- factor(df$ToneDur)

# Run ART separately per tone duration
for(td in levels(df$ToneDur)) {
  cat("\n### Tone Duration:", td, "###\n")
  
  d <- subset(df, ToneDur == td)
  
  m <- art(FTPLrating ~ Predp34 * p34 + Error(Subj/(Predp34*p34)), data=d)
  
  print(anova(m))
}