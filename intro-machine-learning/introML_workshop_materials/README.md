# Introduction to Machine Learning Workshop

## Workshop Description

Machine learning is becoming increasingly important, not only in data science, but also in bioinformatics. Machine learning gives you the opportunity to find hidden patterns in your data and draw predictive conclusions. In this workshop, you'll learn key machine learning concepts and approaches, and get hands-on practice with simple examples.

### Topics Covered

- Supervised and unsupervised machine learning
- Training data and test data sets
- Bias and variance
- Cross validation
- Performance evaluation

### Algorithms

- Random Forest
- K-means clustering

---

## Workshop Materials

This folder contains organized workshop materials:

### Hands-on Exercises

The `hands-on/` folder contains the following R Markdown files:

| File | Description |
|------|-------------|
| `hands-on/01_data_exploration_preprocessing.Rmd` | Data exploration, handling missing values, encoding, scaling |
| `hands-on/02_kmeans.Rmd` | K-means unsupervised clustering |
| `hands-on/03_random_forest.Rmd` | Random Forest classification (supervised learning) |

### Data Files

The `hands-on/data/` folder contains:
- `covid-symptoms.csv` - Original COVID-19 symptoms dataset (input for Part 1)

**Note:** Part 1 (`01_data_exploration_preprocessing.Rmd`) creates `covid_ml_clean.csv` as output, which is then used in Parts 2 and 3.

### Slides

The `slides/` folder contains:
- `IntroML_workshop_011526.pptx` - Workshop presentation slides

---

## Dataset: COVID-19 Symptoms and Test Results

### Description

These workshop materials use a COVID-19 clinical dataset to introduce machine learning concepts. The dataset contains patient symptoms, vitals, and COVID-19 test outcomes.

We use predictor variables (symptoms, comorbidities, vitals) to estimate the probability of an outcome variable (positive or negative COVID-19 test).

**Important note:** Our focus is on building technical skills in machine learning. Although some symptoms or combinations of symptoms may correlate with test results, our goal is not to build a highly predictive model for COVID-19. Instead, we use this dataset to motivate exploration of ML techniques and concepts.

### Features

The dataset includes the following variables:

| Variable | Description |
|----------|-------------|
| `covid19_test_results` | Target variable (Positive/Negative) |
| `temperature` | Patient temperature |
| `high_risk_exposure_occupation` | Binary indicator |
| `high_risk_interactions` | Binary indicator |
| `cough` | Binary indicator |
| `cough_severity` | Categorical (Mild/Moderate/Severe) |
| `fever` | Binary indicator |
| `sob` | Shortness of breath |
| `fatigue` | Binary indicator |
| `headache` | Binary indicator |
| `loss_of_smell` | Binary indicator |
| `loss_of_taste` | Binary indicator |
| `sore_throat` | Binary indicator |
| ... | Additional symptom indicators |

### Data Citation

```
@dataset{2020covidclinicaldata,
  author =       {Carbon Health and Braid Health},
  title =        {Coronavirus Disease 2019 (COVID-19) Clinical Data Repository},
  howpublished = {Accessed from \url{https://covidclinicaldata.org/}},
  year =         2020,
  version =      {10-20-2020}
}
```

---

## Acknowledgments

We thank the **UCSF Library Data Science Initiative** for sharing this dataset with us. Their instructors have developed excellent materials using this data, and we are grateful for their collaboration and willingness to share resources with the bioinformatics community.

Original dataset repository: [Covid-Test-Predictions](https://github.com/geoffswc/Covid-Test-Predictions)

---

## Prerequisites

### Required R Packages

```r
install.packages(c(
  "tidyverse",
  "ggplot2", 
  "caret",
  "randomForest",
  "cluster",
  "factoextra",
  "skimr",
  "naniar"
))
```

### Recommended Background

- Basic familiarity with R programming
- Understanding of basic statistics (mean, standard deviation)
- No prior machine learning experience required!

---

## Contact

Gladstone Bioinformatics Core

---

*Last updated: January 2026*
