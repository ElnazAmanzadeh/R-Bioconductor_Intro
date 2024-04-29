
# Data Handling
# Reading File
#solution1
download.file(url = "https://github.com/carpentries-incubator/bioc-intro/raw/main/episodes/data/rnaseq.csv",
              destfile = "data/rnaseq.csv")
#solution2
setwd("E:/courses/R Workshop/data")
read.table("rnaseq.csv", sep = ",", header=TRUE)

rna <- read.csv("rnaseq.csv")
rna

head(rna)


# DataFrame
# Creating a dataframe
df <- data.frame(
  Name = c("Alice", "Bob", "Charlie", "David"),
  Age = c(25, 30, 35, 40),
  Gender = c("Female", "Male", "Male", "Male")
)

# Print the dataframe
print(df)

#Structure of a dataframe using str
str(df)
str(rna)

dim(rna)
nrow(rna)
ncol(rna)
names(rna)
rownames(rna)
summary(rna)

# Indexing and subsetting a dataframe
# first element in the first column of the data frame (as a vector)
rna[1, 1]
# first element in the 6th column (as a vector)
rna[1, 6]
# first column of the data frame (as a vector)
rna[, 1]
# first column of the data frame (as a data.frame)
rna[1]
# first three elements in the 7th column (as a vector)
rna[1:3, 7]
# the 3rd row of the data frame (as a data.frame)
rna[3, ]
# equivalent to head_rna <- head(rna)
head_rna <- rna[1:6, ]
head_rna


#Exclude
rna[, -1]          ## The whole data frame, except the first column
rna[-c(7:66465), ] ## Equivalent to head(rna)



#Data frames can be subsetted by calling indices 
rna["gene"]       # Result is a data.frame
rna[, "gene"]     # Result is a vector
rna[["gene"]]     # Result is a vector
rna$gene



#Factors represent categorical data

sex <- factor(c("male", "female", "female", "male", "female"))

levels(sex)
nlevels(sex)


sex ## current order


sex <- factor(sex, levels = c("male", "female"))
sex ## after re-ordering

plot(sex)

as.character(sex)
levels(sex)
levels(sex) <- c("M", "F")
levels(sex)
sex

#Creating a dataframe using data.frame() function
animal_data <- data.frame(
  animal = c("dog", "cat", "sea cucumber", "sea urchin"),
  feel = c("furry", "squishy", "spiny", "fluffy"),
  weight = c(45, 8, 1.1, 0.8))

#mistakes?


#Matrices
m <- matrix(1:9, ncol = 3, nrow = 3)
m

#Date: YEAR-MONTH-DAY
library("lubridate")

my_date <- ymd("2015-01-01")
str(my_date)

# sep indicates the character to use to separate each component
my_date <- ymd(paste("2015", "1", "1", sep = "-"))
str(my_date)


x <- data.frame(year = c(1996, 1992, 1987, 1986, 2000, 1990, 2002, 1994, 1997, 1985),
                month = c(2,  3,  3, 10,  1,  8,  3,  4,  5,  5),
                day = c(24,  8,  1,  5,  8, 17, 13, 10, 11, 24),
                value = c(4,  5,  1,  9,  3,  8, 10,  2,  6,  7))
x



paste(x$year, x$month, x$day, sep = "-")


ymd(paste(x$year, x$month, x$day, sep = "-"))



#Lists
l <- list(1:10, ## numeric
          letters, ## character
          installed.packages(), ## a matrix
          cars, ## a data.frame
          list(1, 2, 3)) ## a list
length(l)

l
str(l)

l[[1]] ## first element



l[1:2] ## a list of length 2


l[1]   ## a list of length 1

#Exporting and writing tabular data
write.csv(rna, file = "data_output/my_rna.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Manipulate Data 

#tidyr
#dplyr
#tidyverse is an umbrella-package
## load the tidyverse packages, incl. dplyr
library("tidyverse")

#BiocManager::install("tidyverse")
#Reading Data
rna <- read_csv("rnaseq.csv")

rna

#Selecting columns

select(rna, age, sample, tissue, expression)

#select and unselect
select(rna, -tissue, -organism)

#filter
filter(rna, sex == "Male")

filter(rna, sex == "Male" & infection == "NonInfected")

genes <- select(rna, gene, hsapiens_homolog_associated_gene_name)
genes

filter(genes, is.na(hsapiens_homolog_associated_gene_name))

filter(genes, !is.na(hsapiens_homolog_associated_gene_name))


#Pipes


rna2 <- filter(rna, sex == "Male")
rna3 <- select(rna2, gene, sample, tissue, expression)
rna3


rna3 <- select(filter(rna, sex == "Male"), gene, sample, tissue, expression)
rna3


#ctrl+shift+M = %>% 
rna %>%
  filter(sex == "Male") %>%
  select(gene, sample, tissue, expression)



rna3 <- rna %>%
  filter(sex == "Male") %>%
  select(gene, sample, tissue, expression)

rna3














#Mutate
rna %>%
  mutate(time_hours = time * 24) %>%
  select(time, time_hours)


rna %>%
  mutate(time_hours = time * 24,
         time_mn = time_hours * 60) %>%
  select(time, time_hours, time_mn)






#Split-apply-combine data analysis

rna %>%
  group_by(gene)
rna %>%
  group_by(sample)


rna %>%
  group_by(gene) %>%
  summarise(mean_expression = mean(expression))


rna %>%
  group_by(sample) %>%
  summarise(mean_expression = mean(expression))


rna %>%
  group_by(gene, infection, time) %>%
  summarise(mean_expression = mean(expression),
            median_expression = median(expression))






#Counting

rna %>%
  count(infection)

rna %>%
  group_by(infection) %>%
  summarise(n = n())

rna %>%
  count(infection, time)

rna %>%
  group_by(infection, time) %>%
  summarise(n = n())


rna %>%
  count(infection, time) %>%
  arrange(time)


rna %>%
  count(infection, time) %>%
  arrange(n)


rna %>%
  count(infection, time) %>%
  arrange(desc(n))


#pivot_wider
rna_wide <- rna_exp %>%
  pivot_wider(names_from = sample,
              values_from = expression)
rna_wide







#Reshaping data 


rna %>%
  arrange(gene)


rna_exp <- rna %>%
  select(gene, sample, expression)
rna_exp


#pivot_wider adds NA for missing values

rna_with_missing_values <- rna %>%
  select(gene, sample, expression) %>%
  filter(gene %in% c("Asl", "Apod", "Cyp2d22")) %>%
  filter(sample %in% c("GSM2545336", "GSM2545337", "GSM2545338")) %>%
  arrange(sample) %>%
  filter(!(gene == "Cyp2d22" & sample != "GSM2545338"))
rna_with_missing_values

#parameterising with values_fill in pivot_sider
rna_with_missing_values %>%
  pivot_wider(names_from = sample,
              values_from = expression)

#filling vlues
rna_with_missing_values %>%
  pivot_wider(names_from = sample,
              values_from = expression,
              values_fill = 0)



#pivot_longer uses column names and turins them into new values

rna_long <- rna_wide %>%
  pivot_longer(names_to = "sample",
               values_to = "expression",
               -gene)
rna_long

#choose specific columns
rna_wide %>%
  pivot_longer(names_to = "sample",
               values_to = "expression",
               cols = starts_with("GSM"))


rna_wide %>%
  pivot_longer(names_to = "sample",
               values_to = "expression",
               GSM2545336:GSM2545380)




wide_with_NA <- rna_with_missing_values %>%
  pivot_wider(names_from = sample,
              values_from = expression)
wide_with_NA



wide_with_NA %>%
  pivot_longer(names_to = "sample",
               values_to = "expression",
               -gene)




# Joining Tables
rna_mini <- rna %>%
  select(gene, sample, expression) %>%
  head(10)
rna_mini





















###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data Visualization

# ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) +  <GEOM_FUNCTION>()

ggplot(data = rna)


ggplot(data = rna, mapping = aes(x = expression))


ggplot(data = rna, mapping = aes(x = expression)) +
  geom_histogram()


# Assign plot to a variable
rna_plot <- ggplot(data = rna,
                   mapping = aes(x = expression))

# Draw the plot
rna_plot + geom_histogram()



rna <- rna %>%
  mutate(expression_log = log2(expression + 1))




ggplot(rna, aes(x = expression_log)) + geom_histogram()






# This is the correct syntax for adding layers
rna_plot +
  geom_histogram()

# This will not add the new layer and will return an error message
rna_plot
+ geom_histogram()











# Building your plots iteratively 

rna_fc <- rna %>% select(gene, time,
                         gene_biotype, expression_log) %>%
  group_by(gene, time, gene_biotype) %>%
  summarize(mean_exp = mean(expression_log)) %>%
  pivot_wider(names_from = time,
              values_from = mean_exp) %>%
  mutate(time_8_vs_0 = `8` - `0`, time_4_vs_0 = `4` - `0`)
rna_fc



ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_point()



ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_point(alpha = 0.3)






ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_point(alpha = 0.3, color = "blue")




ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0)) +
  geom_point(alpha = 0.3, aes(color = gene_biotype))






ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0,
                                    color = gene_biotype)) +
  geom_point(alpha = 0.3)




ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0,
                                    color = gene_biotype)) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0)



ggplot(data = rna_fc, mapping = aes(x = time_4_vs_0, y = time_8_vs_0,
                                    color = gene_biotype)) +
  geom_jitter(alpha = 0.3) +
  geom_abline(intercept = 0)














# Boxplot
ggplot(data = rna,
       mapping = aes(y = expression_log, x = sample)) +
  geom_boxplot()



ggplot(data = rna,
       mapping = aes(y = expression_log, x = sample)) +
  geom_jitter(alpha = 0.2, color = "tomato") +
  geom_boxplot(alpha = 0)





ggplot(data = rna,
       mapping = aes(y = expression_log, x = sample)) +
  geom_jitter(alpha = 0.2, color = "tomato") +
  geom_boxplot(alpha = 0) +
  theme(axis.text.x = element_text(angle = 90,  hjust = 0.5, vjust = 0.5))






#Line plots

rna_fc <- rna_fc %>% arrange(desc(time_8_vs_0))

genes_selected <- rna_fc$gene[1:10]
genes_selected

sub_rna <- rna %>%
  filter(gene %in% genes_selected)
sub_rna


mean_exp_by_time <- sub_rna %>%
  group_by(gene,time) %>%
  summarize(mean_exp = mean(expression_log))


mean_exp_by_time



ggplot(data = mean_exp_by_time, mapping = aes(x = time, y = mean_exp)) +
  geom_line()


ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp, group = gene)) +
  geom_line()





ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp, color = gene)) +
  geom_line()














#Faceting
ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp)) + geom_line() +
  facet_wrap(~ gene)




ggplot(data = mean_exp_by_time,
       mapping = aes(x = time, y = mean_exp)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y")




mean_exp_by_time_sex <- sub_rna %>%
  group_by(gene, time, sex) %>%
  summarize(mean_exp = mean(expression_log))




ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y")






ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank())



# One column, facet by rows
ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = gene)) +
  geom_line() +
  facet_grid(sex ~ .)








# One row, facet by column
ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = gene)) +
  geom_line() +
  facet_grid(. ~ sex)
















# Customisation

ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Mean gene expression by duration of the infection",
       x = "Duration of the infection (in days)",
       y = "Mean gene expression")





















ggplot(data = mean_exp_by_time_sex,
       mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Mean gene expression by duration of the infection",
       x = "Duration of the infection (in days)",
       y = "Mean gene expression")  +
  theme(text = element_text(size = 16))





















ggplot(data = mean_exp_by_time_sex,
      mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Mean gene expression by duration of the infection",
       x = "Duration of the infection (in days)",
       y = "Mean gene expression")  +
  theme(text = element_text(size = 16),
        axis.text.x = element_text(colour = "royalblue4", size = 12),
        axis.text.y = element_text(colour = "royalblue4", size = 12),
        panel.grid = element_line(colour="lightsteelblue1"),
        legend.position = "top")






















blue_theme <- theme(axis.text.x = element_text(colour = "royalblue4",
                                               size = 12),
                    axis.text.y = element_text(colour = "royalblue4",
                                               size = 12),
                    text = element_text(size = 16),
                    panel.grid = element_line(colour="lightsteelblue1"))

ggplot(rna, aes(x = expression_log)) +
  geom_histogram(bins = 20) +
  blue_theme












# Composing plots

rna$chromosome_name <- factor(rna$chromosome_name,
                              levels = c(1:19,"X","Y"))

count_gene_chromosome <- rna %>% select(chromosome_name, gene) %>%
  distinct() %>% ggplot() +
  geom_bar(aes(x = chromosome_name), fill = "seagreen",
           position = "dodge", stat = "count") +
  labs(y = "log10(n genes)", x = "chromosome") +
  scale_y_log10()

count_gene_chromosome




















exp_boxplot_sex <- ggplot(rna, aes(y=expression_log, x = as.factor(time),
                                   color=sex)) +
  geom_boxplot(alpha = 0) +
  labs(y = "Mean gene exp",
       x = "time") + theme(legend.position = "none")

exp_boxplot_sex

















#install.packages("patchwork")

library("patchwork")
count_gene_chromosome + exp_boxplot_sex

count_gene_chromosome / exp_boxplot_sex



count_gene_chromosome + exp_boxplot_sex + plot_layout(ncol = 1)




count_gene_chromosome +
  (count_gene_chromosome + exp_boxplot_sex) +
  exp_boxplot_sex +
  plot_layout(ncol = 1)




count_gene_chromosome /
  (count_gene_chromosome | exp_boxplot_sex) /
  exp_boxplot_sex


















install.packages("gridExtra")


library("gridExtra")
grid.arrange(count_gene_chromosome, exp_boxplot_sex, ncol = 2)























# Exporting Plots
my_plot <- ggplot(data = mean_exp_by_time_sex,
                  mapping = aes(x = time, y = mean_exp, color = sex)) +
  geom_line() +
  facet_wrap(~ gene, scales = "free_y") +
  labs(title = "Mean gene expression by duration of the infection",
       x = "Duration of the infection (in days)",
       y = "Mean gene expression") +
  guides(color=guide_legend(title="Gender")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "royalblue4", size = 12),
        axis.text.y = element_text(colour = "royalblue4", size = 12),
        text = element_text(size = 16),
        panel.grid = element_line(colour="lightsteelblue1"),
        legend.position = "top")
ggsave("fig_output/mean_exp_by_time_sex.png", my_plot, width = 15,
       height = 10)

# This also works for grid.arrange() plots
combo_plot <- grid.arrange(count_gene_chromosome, exp_boxplot_sex,
                           ncol = 2, widths = c(4, 6))
ggsave("fig_output/combo_plot_chromosome_sex.png", combo_plot,
       width = 10, dpi = 300)
