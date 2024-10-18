# Statistics overview




# Table of contents
* [Resources](#resources)
* [Performing statistical tests](#performing-statistical-tests)
  * [Types of variables](#types-of-variables)
* [Choosing the correct test](#choosing-the-correct-test)
* [Statistical_tests](#statistical-tests)


---

# Resources

* [GUide to STatistical Analysis in Microbial Ecology (GUSTA ME)](https://mb3is.megx.net/gustame)
* [LHS 610: Exploratory Data Analysis for Health](https://kdpsingh.lab.medicine.umich.edu/lhs-610)
  * We haven't personally tried this course, but they provide great videos and code examples for learning how to explore data using R.
* [#bioinformatics live twitter feed](https://twitter.com/search?q=%23bioinformatics&src=hash)
* [ggpubr](http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/)
  * Nice package for "publication-ready" figures.
* [Harvard's Data Science: R Basics](https://www.edx.org/course/data-science-r-basics)
* [Collaborative spreadsheet of resources](https://docs.google.com/document/d/1A9BbOCsrg1ikLaBltKhXVKj-eetlrBqR-1u-2V99I2c/edit#)
* [Choose the right test](https://stats.idre.ucla.edu/other/mult-pkg/whatstat/)


---

# Performing statistical tests

You can perform statistical tests on data that have been collected in a statistically valid manner – either through an experiment, or through observations made using probability sampling methods.

For a statistical test to be valid, your sample size needs to be large enough to approximate the true distribution of the population being studied.

To determine which statistical test to use, you need to know:
* whether your data meets certain assumptions.
* the types of variables that you’re dealing with.


## Statistical assumptions

Statistical tests make some common assumptions about the data they are testing:

* Independence of observations (a.k.a. no autocorrelation): The observations/variables you include in your test are not related (for example, multiple measurements of a single test subject are not independent, while measurements of multiple different test subjects are independent).
* Homogeneity of variance: the variance within each group being compared is similar among all groups. If one group has much more variation than others, it will limit the test’s effectiveness.
* Normality of data: the data follows a normal distribution (a.k.a. a bell curve). This assumption applies only to quantitative data.

If your data do not meet the assumptions of normality or homogeneity of variance, you may be able to perform a nonparametric statistical test, which allows you to make comparisons without any assumptions about the data distribution.

If your data do not meet the assumption of independence of observations, you may be able to use a test that accounts for structure in your data (repeated-measures tests or tests that include blocking variables).

## Types of variables
The types of variables you have usually determine what type of statistical test you can use.

Quantitative variables represent amounts of things (e.g. the number of trees in a forest). Types of quantitative variables include:

* Continuous (a.k.a ratio variables): represent measures and can usually be divided into units smaller than one (e.g. 0.75 grams).
* Discrete (a.k.a integer variables): represent counts and usually can’t be divided into units smaller than one (e.g. 1 tree).

Categorical variables represent groupings of things (e.g. the different tree species in a forest). Types of categorical variables include:

* Ordinal: represent data with an order (e.g. rankings).
* Nominal: represent group names (e.g. brands or species names).
* Binary: represent data with a yes/no or 1/0 outcome (e.g. win or lose).

Choose the test that fits the types of predictor and outcome variables you have collected (if you are doing an experiment, these are the independent and dependent variables). Consult the tables below to see which test best matches your variables.


--- 

# Choosing the correct test

## Choosing a parametric test: regression, comparison, or correlation

Parametric tests usually have stricter requirements than nonparametric tests, and are able to make stronger inferences from the data. They can only be conducted with data that adheres to the common assumptions of statistical tests.

The most common types of parametric test include regression tests, comparison tests, and correlation tests. Here's a flow chart to help you pick the right test.

![img ](https://cdn.scribbr.com/wp-content/uploads/2020/01/flowchart-for-choosing-a-statistical-test-993x1024.png)


### Regression tests

Regression tests are used to test cause-and-effect relationships. They look for the effect of one or more continuous variables on another variable.

|		|	Predictor variable	|	Outcome variable	|	Research question example	|
|	--------	|	--------	|	--------	|	--------	|
|	Simple linear regression	|	Continuous 1 predictor	|	Continuous 1 outcome	|	What is the effect of income on longevity?	|
|	Multiple linear regression	|	Continuous 2 or more predictors	|	Continuous 1 outcome	|	What is the effect of income and minutes of exercise per day on longevity?	|
|	Logistic regression	|	Continuous	|	Binary	|	What is the effect of drug dosage on the survival of a test subject?	|


### Comparison tests

Comparison tests look for differences among group means. They can be used to test the effect of a categorical variable on the mean value of some other characteristic.

T-tests are used when comparing the means of precisely two groups (e.g. the average heights of men and women). ANOVA and MANOVA tests are used when comparing the means of more than two groups (e.g. the average heights of children, teenagers, and adults).


| Predictor variable |         Outcome variable        |              Research question example              |                                                                                                                     |   |
|:------------------:|:-------------------------------:|:---------------------------------------------------:|---------------------------------------------------------------------------------------------------------------------|---|
| Paired t-test      | Categorical 1 predictor         | Quantitative groups come from the same population   | What is the effect of two different test prep programs on the average exam scores for students from the same class? |   |
| Independent t-test | Categorical 1 predictor         | Quantitative groups come from different populations | What is the difference in average exam scores for students from two different schools?                              |   |
| ANOVA              | Categorical 1 or more predictor | Quantitative 1 outcome                              | What is the difference in average pain levels among post-surgical patients given three different painkillers?       |   |
| MANOVA             | Categorical 1 or more predictor | Quantitative 2 or more outcome                      | What is the effect of flower species on petal length, petal width, and stem length?                                 |   |



### Correlation tests

Correlation tests check whether two variables are related without assuming cause-and-effect relationships.

These can be used to test whether two variables you want to use in (for example) a multiple regression test are autocorrelated.

|            | Predictor variable | Outcome variable |                                      Research question example                                     |
|------------|:------------------:|:----------------:|:--------------------------------------------------------------------------------------------------:|
| Pearson    |     Continuous     |    Continuous    | How are latitude and temperature related?                                                          |
| Chi-Square |     Categorical    |    Categorical   | How is membership in a sports team related to membership in drama club among high school students? |


## Choosing a nonparametric test

Non-parametric tests don’t make as many assumptions about the data, and are useful when one or more of the common statistical assumptions are violated. However, the inferences they make aren’t as strong as with parametric tests.

|                           |      Predictor variable      |                   Outcome variable                  |         Use in place of…         |
|---------------------------|:----------------------------:|:---------------------------------------------------:|:--------------------------------:|
| Spearman                  | Ordinal                      | Ordinal                                             | Regression and correlation tests |
| Sign test                 | Categorical                  | Quantitative                                        | T-test                           |
| Kruskal–Wallis            | Categorical 3 or more groups | Quantitative                                        | ANOVA                            |
| ANOSIM                    | Categorical 3 or more groups | Quantitative 2 or more outcome variables            | MANOVA                           |
| Wilcoxon Rank-Sum test    | Categorical 2 groups         | Quantitative groups come from different populations | Independent t-test               |
| Wilcoxon Signed-rank test | Categorical 2 groups         | Quantitative groups come from the same population   | Paired t-test                    |





---

## Statistical tests

## One sample t-test
A one sample t-test allows us to test whether a sample mean (of a normally distributed interval variable) significantly differs from a hypothesized value. 

## One sample median test
Wilcoxon test
A one sample median test allows us to test whether a sample median differs significantly from a hypothesized value.  We will use the same variable, write, as we did in the one sample t-test example above, but we do not need to assume that it is interval and normally distributed (we only need to assume that write is an ordinal variable).

## Binomial test
A one sample binomial test allows us to test whether the proportion of successes on a two-level categorical dependent variable significantly differs from a hypothesized value.

## Chi-square goodness of fit
A chi-square goodness of fit test allows us to test whether the observed proportions for a categorical variable differ from hypothesized proportions. 


## Two independent samples t-test
An independent samples t-test is used when you want to compare the means of a normally distributed interval dependent variable for two independent groups. Use "equal variances" if applicable.


## Wilcoxon-Mann-Whitney test
The Wilcoxon-Mann-Whitney test is a non-parametric analog to the independent samples t-test and can be used when you do not assume that the dependent variable is a normally distributed interval variable (you only assume that the variable is at least ordinal).

## Chi-square test
A chi-square test is used when you want to see if there is a relationship between two categorical variables.

## Fisher’s exact test
The Fisher’s exact test is used when you want to conduct a chi-square test but one or more of your cells has an expected frequency of five or less.  Remember that the chi-square test assumes that each cell has an expected frequency of five or more, but the Fisher’s exact test has no such assumption and can be used regardless of how small the expected frequency is. 

## One-way ANOVA
A one-way analysis of variance (ANOVA) is used when you have a categorical independent variable (with two or more categories) and a normally distributed interval dependent variable and you wish to test for differences in the means of the dependent variable broken down by the levels of the independent variable.

## Kruskal Wallis test
The Kruskal Wallis test is used when you have one independent variable with two or more levels and an ordinal dependent variable. In other words, it is the non-parametric version of ANOVA and a generalized form of the Mann-Whitney test method since it permits two or more groups. If some of the scores receive tied ranks, then a correction factor is used, yielding a slightly different value of chi-squared. 

## Paired t-test
A paired (samples) t-test is used when you have two related observations (i.e., two observations per subject) and you want to see if the means on these two normally distributed interval variables differ from one another. 

## Wilcoxon signed rank sum test
The Wilcoxon signed rank sum test is the non-parametric version of a paired samples t-test.  You use the Wilcoxon signed rank sum test when you do not wish to assume that the difference between the two variables is interval and normally distributed (but you do assume the difference is ordinal). 

## McNemar test
You would perform McNemar’s test if you were interested in the marginal frequencies of two binary outcomes. These binary outcomes may be the same outcome variable on matched pairs (like a case-control study) or two outcome variables from a single group.

## One-way repeated measures ANOVA
You would perform a one-way repeated measures analysis of variance if you had one categorical independent variable and a normally distributed interval dependent variable that was repeated at least twice for each subject.  This is the equivalent of the paired samples t-test, but allows for two or more levels of the categorical variable. This tests whether the mean of the dependent variable differs by the categorical variable. 

## Repeated measures logistic regression
If you have a binary outcome measured repeatedly for each subject and you wish to run a logistic regression that accounts for the effect of multiple measures from single subjects, you can perform a repeated measures logistic regression.


---

