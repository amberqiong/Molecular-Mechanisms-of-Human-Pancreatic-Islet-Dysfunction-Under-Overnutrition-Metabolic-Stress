INS=read_csv("/Users/liguo/Desktop/lipoglucotoxicity/snow/INS.csv")

rawdata<- INS %>%
  filter(condition %in% c("3.3mM","16.7mM"))

rawdata <- rawdata %>%
  mutate(treatment_condition=interaction(treatment,condition))

summary_data <- rawdata %>%
  group_by(treatment_condition) %>%
  summarise(mean_INS = mean(INS), 
            se_INS = sd(INS) / sqrt(n()))

summary_data$treatment_condition <- factor(summary_data$treatment_condition,
                                           levels = c("Ctrl.3.3mM", "Ctrl.16.7mM", "GL.3.3mM", "GL.16.7mM"))



p1=ggplot(summary_data, aes(x = treatment_condition, y = mean_INS,fill=treatment_condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.6, color="black",size=0.8) +
  scale_fill_manual(values = c("#d9f0e8", "#7dccba", "#dbd9f2", "#a39cc2"))+
  geom_errorbar(aes(ymin = mean_INS - se_INS, ymax = mean_INS + se_INS), 
                width = 0.2, size=1, position = position_dodge(width = 0.9)) +
  geom_jitter(data = rawdata, aes(x = treatment_condition, y = INS),
              width = 0.2, shape = 16, size = 2, color = "black")+ylim(0,60)+
  theme_classic()


rawdata2<- INS %>%
  filter(condition %in% c("S.I."))

summary_data2 <- rawdata2 %>%
  group_by(treatment) %>%
  summarise(mean_INS = mean(INS), 
            se_INS = sd(INS) / sqrt(n()))

summary_data2$treatment <- factor(summary_data2$treatment,
                                  levels = c("Ctrl","GL"))



p2=ggplot(summary_data2, aes(x = treatment, y = mean_INS,fill=treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.6, color="black") +
  scale_fill_manual(values = c("#7dccba","#a39cc2"))+
  geom_errorbar(aes(ymin = mean_INS - se_INS, ymax = mean_INS + se_INS), 
                width = 0.2,position = position_dodge(width = 0.9)) +
  geom_jitter(data = rawdata2, aes(x = treatment, y = INS),
              width = 0.2, shape = 16, size = 2.5, color = "black")+
  theme_classic()+ylim(0,4)


viability=read.csv("/Users/liguo/Desktop/lipoglucotoxicity/snow/viability.csv")

viability_2 <- viability %>%
  group_by(treatment) %>%
  summarise(mean_viability = mean(viability), 
            se_viability = sd(viability) / sqrt(n()))

viability_2$treatment <- factor(viability_2$treatment,
                                  levels = c("Ctrl","GL"))

p3=ggplot(viability_2, aes(x = treatment, y = mean_viability,fill=treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),width = 0.6, color="black") +
  scale_fill_manual(values = c("#7dccba","#a39cc2"))+
  geom_errorbar(aes(ymin = mean_viability - se_viability, ymax = mean_viability + se_viability), 
                width = 0.2,position = position_dodge(width = 0.9)) +
  geom_jitter(data = viability, aes(x = treatment, y = viability),
              width = 0.2, shape = 16, size = 2.5, color = "black")+
  theme_classic()+ylim(0,100)
