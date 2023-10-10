
#####PACKAGES
# install.packages("tidyverse")
library(tidyverse) 
# install.packages("dplyr")
library(dplyr)
# install.packages("ggplot2")
library(ggplot2) 

#####DOSAGE RESPONSE -----

#Organize Data
df2p <- read.csv("CO2.csv", header = , stringsAsFactors = FALSE)
colnames(df2p)[1] <- "P.CO2"
outlier <- df2p %>% group_by(P.CO2, Blocker) %>% 
  dplyr::summarize(out.high=mean(Exudate)+2.5*sd(Exudate), out.low=mean(Exudate)-2.5*sd(Exudate))
df2a<- left_join(df2p, outlier, by= c("P.CO2", "Blocker"))
df2 <- df2a[!df2a$Exudate>df2a$out.high | df2a$Exudate<df2a$out.low, ]

#Summarize Data
sum.stats2 <- df2 %>% 
  group_by(P.CO2, Blocker) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(Exudate, na.rm = TRUE),
    sd = sd(Exudate, na.rm = TRUE),
    sem = sd/sqrt(n) )

#Graph: Just CO2
sum.stats2 %>% 
  filter(Blocker == "None") %>%
  ggplot( aes(x=P.CO2, y=mean, color=Blocker, group=Blocker)) + 
  geom_line(size = .5, alpha=1) +
  geom_point(size = 2,alpha=1 ) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),width=.5,lwd=.7,alpha=1) +
  scale_color_manual(values = c("#443F89")) +
  labs(x = "% CO2", y = "Exudate (mg)") +
  theme_classic() +
  ylim(5, 150)

#Graph: CO2 and Water
sum.stats2 %>% 
  ggplot( aes(x=P.CO2, y=mean, color=Blocker, group=Blocker)) + 
  geom_line(aes(linetype=Blocker),size = .75, alpha=1) +
  geom_point(size = 3,alpha=1 ) +
  scale_color_manual(values = c("#70C678","#008791"), name="Treatment", labels=c("No Treatment", "Water")) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),width=.5,lwd=.7,alpha=1, size=.75) +
  labs(x = "% CO2", y = "Exudate (mg)") +
  theme_classic() +
  theme(text=element_text(size=18),
        axis.text.x = element_text(size=16)) +
  ylim(5, 150)

ggsave("Figure1.pdf",units = "in",height = 9,width = 9)

#Statistics
df2$P.CO2 <- as.factor(df2$P.CO2)
df2$Blocker <- as.factor(df2$Blocker)
summary(aov(Exudate ~ P.CO2+Blocker , data = df2))
TukeyHSD(aov(Exudate ~ P.CO2*Blocker, data = df2))

#####BLOCKERS ----

#Organize Data
df <- read.csv("blockers.csv", header = , stringsAsFactors = FALSE)
colnames(df)[1] <- "P.CO2"
df$P.CO2 <- factor(df$P.CO2, levels = c("0", "50", "100"))
df$Blocker <- factor(df$Blocker, levels = c("Ruthenium Red", "Amiloride", "ZnCl2", "HC30013", "Diminazene Aceturate", "Methylene Blue", "Acetazolamide", "Indisulam", "U-104", "S4 ", "Topiramate"))
df$Vehicle <- factor(df$Vehicle, levels = c("Water", "Methyl Cellulose", "1% DMSO"))

#Summarize Data
sum.stats <- df %>% 
  group_by(P.CO2, Blocker, Concentration, Treatment, Vehicle) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(Exudate, na.rm = TRUE),
    sd = sd(Exudate, na.rm = TRUE),
    sem = sd/sqrt(n) )

#Carbonic Anhydrase Inhibitor Graphs
sum.stats %>% 
  filter(Blocker == "Acetazolamide" | Blocker == "Indisulam" | Blocker == "S4 " |
           Blocker == "Topiramate" | Blocker == "U-104") %>%
  ggplot(aes(x=P.CO2, y=mean, color = Treatment, shape = Vehicle, group=interaction(Blocker, Concentration, Treatment, Vehicle))) + 
  geom_line(aes(linetype=Treatment),linewidth = .75, alpha=1) +
  geom_point(size = 3,alpha=1) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, color = Treatment),width=.5,lwd=.7,alpha=1, size=.75) +
  scale_color_manual(values = c("#70C678","#008791")) +
  ylim(-5, 200) +
  facet_wrap(vars(Blocker)) +
  labs(x = "% CO2", y = "Exudate (mg)") +
  theme_classic() +
  theme(text=element_text(size=18),
        axis.text.x = element_text(size=16))

ggsave("Figure7A.pdf",units = "in",height = 9,width = 11)

#Other Blocker Graphs
sum.stats %>% 
  filter(Blocker == "Ruthenium Red" | Blocker == "HC30013" | Blocker == "Amiloride" |
           Blocker == "Diminazene Aceturate" | Blocker == "Methylene Blue" | Blocker == "ZnCl2") %>%
  ggplot(aes(x=P.CO2, y=mean, color = Treatment, shape = Vehicle, group=interaction(Blocker, Concentration, Treatment, Vehicle))) + 
  geom_line(aes(linetype=Treatment),linewidth = .75, alpha=1) +
  geom_point(size = 3,alpha=1) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem, color = Treatment),width=.5,lwd=.7,alpha=1, size=.75) +
  scale_color_manual(
    values = c("#70C678", "#97a13b", "#008791"),
  ) +
  ylim(-5, 200) +
  facet_wrap(vars(Blocker)) +
  labs(x = "% CO2", y = "Exudate (mg)") +
  theme_classic() +
  theme(text=element_text(size=18),
        axis.text.x = element_text(size=16)) 

ggsave("Figure7B.pdf",units = "in",height = 9,width = 11)

#Statistics
df.AZ <- df[df$Blocker == "Acetazolamide",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.AZ))
df.AZ$P.CO2 <- as.factor(df.AZ$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.AZ))

df.AM <- df[df$Blocker == "Amiloride",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.AM))
df.AM$P.CO2 <- as.factor(df.AM$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.AM))

df.DA <- df[df$Blocker == "Diminazene Aceturate",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.DA))
df.DA$P.CO2 <- as.factor(df.DA$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.DA))

df.HC <- df[df$Blocker == "HC30013",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.HC))
df.HC$P.CO2 <- as.factor(df.HC$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.HC))

df.Ind <- df[df$Blocker == "Indisulam",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.Ind))
df.Ind$P.CO2 <- as.factor(df.Ind$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.Ind))

df.MB <- df[df$Blocker == "Methylene Blue",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.MB))
df.MB$P.CO2 <- as.factor(df.MB$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.MB))

df.RuR <- df[df$Blocker == "Ruthenium Red",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.RuR))
df.RuR$P.CO2 <- as.factor(df.RuR$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.RuR))

df.S4 <- df[df$Blocker == "S4 ",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.S4))
df.S4$P.CO2 <- as.factor(df.S4$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.S4))

df.Top <- df[df$Blocker == "Topiramate",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.Top))
df.Top$P.CO2 <- as.factor(df.Top$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.Top))

df.U104 <- df[df$Blocker == "U-104",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.U104))
df.U104$P.CO2 <- as.factor(df.U104$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.U104))

df.Zn <- df[df$Blocker == "ZnCl2",]
summary(aov(Exudate ~ P.CO2+Concentration , data = df.Zn))
df.Zn$P.CO2 <- as.factor(df.Zn$P.CO2)
TukeyHSD(aov(Exudate ~ P.CO2*Concentration , data = df.Zn))

# Note: For the above data, "0mM" of each blocker is the respective vehicle control. Blockers with the same vehicle are compared to the same control data.

# Bonferroni correction was made "by hand" by dividing alpha by the number of comparisons in each subset of tests

#####AITC CONTROLS ----

#Organize Data
df3 <- read.csv("AITC.csv", header = , stringsAsFactors = FALSE)
colnames(df3)[1] <- "Stimulus"
unique(df3$Treatment)
df3$Treatment <- factor(df3$Treatment, 
  levels = c("Control",  "H2O (Vehicle)", "HC030013", "Acetazolamide", 
    "Amiloride (1 mM)", "Amiloride (5 mM)", 
    "Diminazene Aceturate (0.1 mM)", "Diminazene Aceturate (0.05 mM)","Indisulam")
)
unique(df3$Stimulus)
df3$Stimulus <- factor(df3$Stimulus, levels = c("Mineral Oil (Vehicle)","AITC"))

#Summarize Data
sum.stats <- df3 %>% 
  group_by(Stimulus,Treatment) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(Exudate, na.rm = TRUE),
    sd = sd(Exudate, na.rm = TRUE),
    sem = sd/sqrt(n) )

#AITC Graph
sum.stats %>% 
 ggplot(aes(x=Treatment, y=mean, fill=Stimulus)) +
  geom_bar(stat="identity",  
           position = position_dodge2(width=.75, preserve = "single")) + 
  theme_classic() + 
  xlab("Treatment") + ylab("Exudate (mg)") +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
            lwd=.7,alpha=1,
            position = position_dodge2(width=.75,  preserve = "single")) +
  scale_fill_manual(values = c("#008791","#97a13b"))+
  xlab("")+
  theme(text=element_text(size=18),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=16))

ggsave("Figure9.pdf",units = "in",height = 9,width = 11)

#Statistics
summary(aov(Exudate ~ Stimulus * Treatment , data = df3))
TukeyHSD(aov(Exudate ~ Treatment, data = df3))

df3.stats <- rstatix::tukey_hsd(aov(Exudate ~ Stimulus * Treatment, data = df3))

##### ORGANIC ACIDS ----

df <- read.csv("OrganicAcids.csv")
colnames(df)[1] <- "ChemicalName"

df$ChemicalName <- factor(df$ChemicalName,
                          levels = c("Water Control","Formic acid","Acetic acid","Propionic acid"))
df$Blocker <- factor(df$Blocker,levels=c("Vehicle","Amiloride"))
df$Blocker.concentration <- factor(df$blocker.concentration,c("0 mM","1 mM","5 mM","10 mM"))

# summarize data ---

summary <- df %>% 
  group_by(ChemicalName, Blocker, Blocker.concentration, PPB) %>% 
  summarise(
    mean = mean(Exudate.mg, na.rm = T),
    sd = sd(Exudate.mg, na.rm = T),
    sem = sd(Exudate.mg, na.rm = T)/sqrt(n()),
    n=n()
  )

colnames(summary)[5] <- "Exudate.mg"

# plots ---

summary[summary$Blocker == "Vehicle" & summary$PPB != 0,] %>%
  ggplot( aes(PPB, Exudate.mg, shape = ChemicalName, color = ChemicalName)) +
  geom_line( lwd = 1.5, alpha=1 ) +
  geom_errorbar(aes(ymin = Exudate.mg-sem, ymax = Exudate.mg+sem), width=0.2,color = "black",lwd=1.25) +
  geom_point(position = position_dodge(0.8), size=+4) + 
  labs( x = "Concentration (PPB)", y = "Exudate (mg)") +
  scale_color_manual(values=c( "#7a5195", "#ef5675","#003f5c")) +
  scale_x_log10(limits = c(12,1100)) + annotation_logticks(sides = "b") +
  facet_wrap( ChemicalName ~ .) +
  theme_classic() +
  theme(text = element_text(size = 18),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=16))

ggsave("Figure 8A - OrganicAcids.pdf",units = "in",width = 9,height = 4)

summary[summary$PPB == 0|summary$PPB == 250, ] %>% 
  ggplot( aes( Blocker.concentration, Exudate.mg, fill = ChemicalName),width=0.9) +
  geom_col(position = position_dodge(0.8, preserve = "single"), width = 0.8) +
  geom_errorbar( aes(ymin = Exudate.mg-sem, ymax = Exudate.mg+sem),
                 width = 0.3, position = position_dodge(0.8, preserve = "single"), lwd=+1 ) +
  labs( x = "Amiloride Treatment", y = "Exudate (mg)") +
  scale_fill_manual("Stimulus",values=c("dodgerblue" ,"#7a5195", "#ef5675","#003f5c")) +
  scale_y_continuous(limits = c(0,300)) +
  # facet_wrap( ChemicalName ~ ., nrow = 1) +
  theme_classic( base_size = 10)

ggsave("Figure 8B - OrganicAcids-Amiloride.pdf",units = "in",width = 9,height = 4)

### stats
rstatix::anova_test(aov(Exudate.mg ~ ChemicalName * Blocker.concentration, data = df))
rstatix::tukey_hsd(aov(data = df[df$blocker.concentration == "0 mM",], Exudate.mg ~ ChemicalName))
rstatix::tukey_hsd(aov(data = df[df$blocker.concentration == "1 mM",], Exudate.mg ~ ChemicalName))
rstatix::tukey_hsd(aov(data = df[df$blocker.concentration == "5 mM",], Exudate.mg ~ ChemicalName))
rstatix::tukey_hsd(aov(data = df[df$blocker.concentration == "10 mM",], Exudate.mg ~ ChemicalName))

