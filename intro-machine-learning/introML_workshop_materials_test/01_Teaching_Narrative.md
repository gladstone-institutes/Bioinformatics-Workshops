# Teaching Narrative: Part 1 - Data Exploration & Preprocessing

## Introduction (5 minutes)

**Opening Statement:**
"Welcome to Part 1 of our Introduction to Machine Learning workshop. Today, we're going to learn about one of the most critical but often overlooked aspects of machine learning: data preprocessing. Before we can build any ML model, we need to understand and prepare our data properly."

**Key Points to Emphasize:**
- Data preprocessing is often 80% of the work in ML projects
- Different algorithms have different requirements
- Garbage in = garbage out - poor data quality leads to poor models

---

## Section 1: Load Required Libraries (2 minutes)

**What to Say:**
"Let's start by loading the libraries we'll need. We're using `tidyverse` which includes `dplyr` for data manipulation and `ggplot2` for visualization, and `skimr` for comprehensive data summaries."

**Explain the Code:**
```r
library(tidyverse)  # Loads multiple packages including dplyr, ggplot2, readr, etc.
library(skimr)      # Provides the skim() function for detailed data summaries
```

**Teaching Points:**
- `tidyverse` is a collection of R packages designed for data science
- `skimr` provides better summaries than base R's `summary()` function
- Always load libraries at the beginning of your script

---

## Section 2: Load and Explore the Data (10 minutes)

**What to Say:**
"We're working with a COVID-19 symptoms dataset. This contains patient symptoms and test results - perfect for learning ML because it has both features (symptoms) and a target variable (test results)."

### Step 1: Load the Data
```r
covid_data <- read.csv("data/covid-symptoms.csv")
```

**Explain:**
- `read.csv()` reads a CSV file into R as a data frame
- The file path is relative to your working directory
- The data is stored in the variable `covid_data`

### Step 2: Check Target Distribution
```r
print(table(covid_data$covid19_test_results))
```

**Explain:**
- `table()` creates a frequency table showing counts of each category
- The `$` operator accesses a column in a data frame
- This shows us the distribution of our target variable (what we want to predict)
- **Important:** Check for class imbalance - if one class dominates, we may need special techniques

### Step 3: Comprehensive Data Overview
```r
skim(covid_data)
```

**Explain the `skim()` Function:**
- **Purpose:** Provides a comprehensive overview of your dataset
- **What it shows:**
  - **Dimensions:** Number of rows and columns
  - **Column types:** numeric, character, factor, etc.
  - **Summary statistics:** mean, median, quartiles for numeric columns
  - **Missing values:** `n_missing` shows count, `complete_rate` shows proportion (1.0 = no missing)
  - **Character/factor summaries:** unique values, top values

**Key Observations to Highlight:**
- Point out the `n_missing` column - shows which variables have missing data
- `complete_rate` of 1.0 means no missing values; lower values indicate data quality issues
- Identify which columns are the target (what we predict) vs features (what we use to predict)

---

## Section 3: Handling Missing Data (15 minutes)

**What to Say:**
"Most machine learning algorithms cannot handle missing values. We have two main strategies: remove rows with missing data, or impute (fill in) the missing values. Today we'll use imputation because we don't want to lose data."

### Why Missing Data Matters
- Most algorithms will fail or produce errors with missing values
- Removing rows can significantly reduce dataset size
- Imputation preserves data but requires careful consideration

### Step 1: Create a Copy
```r
covid_imputed <- covid_data
```

**Explain:**
- Always work on a copy to preserve original data
- Good practice for data analysis

### Step 2: Identify Numeric Columns
```r
numeric_cols <- names(covid_imputed)[sapply(covid_imputed, is.numeric)]
```

**Break Down This Line:**
- `is.numeric()` checks if each column is numeric
- `sapply()` applies `is.numeric()` to each column, returns a logical vector
- `names(covid_imputed)[...]` extracts column names where the condition is TRUE
- Result: vector of numeric column names

**Why Exclude Target?**
- We're only imputing features, not the target variable
- Missing targets would mean we can't use those rows for supervised learning anyway

### Step 3: Median Imputation Loop
```r
for (col in numeric_cols) {
  if (sum(is.na(covid_imputed[[col]])) > 0) {
    median_val <- median(covid_imputed[[col]], na.rm = TRUE)
    covid_imputed[[col]][is.na(covid_imputed[[col]])] <- median_val
    cat("Imputed", col, "with median:", median_val, "\n")
  }
}
```

**Explain Each Part:**

**The `for` loop:**
- Iterates through each numeric column
- `col` is the column name (string)

**`is.na()` function:**
- Returns TRUE/FALSE for each value indicating if it's missing
- `sum(is.na(...))` counts how many missing values exist
- `> 0` checks if there are any missing values

**`median()` function:**
- **Parameters:**
  - `x`: The vector/column to calculate median for
  - `na.rm = TRUE`: Remove NA values before calculating (required!)
- **Why median over mean?** Median is robust to outliers

**Double bracket notation `[[col]]`:**
- `covid_imputed[[col]]` accesses a column by name stored in a variable
- Alternative to `covid_imputed$column_name` when the name is in a variable

**Assignment with logical indexing:**
- `covid_imputed[[col]][is.na(...)]` selects only missing values
- Replaces them with the median value

**`cat()` function:**
- Prints messages to console
- Useful for tracking what's happening

### Step 4: Verify No Missing Values Remain
```r
cat("\nRemaining missing values:", sum(is.na(covid_imputed)), "\n")
```

**Explain:**
- `sum(is.na(covid_imputed))` counts ALL missing values in entire dataframe
- Should be 0 after imputation
- Always verify your preprocessing steps!

**Alternative Imputation Methods to Mention:**
- Mean imputation (sensitive to outliers)
- Mode imputation (for categorical)
- K-nearest neighbors imputation (more sophisticated)
- Forward/backward fill (for time series)

---

## Section 4: Handling Categorical Variables (10 minutes)

**What to Say:**
"Some algorithms require numeric input. We need to convert categorical variables to numbers. There are different encoding methods - today we'll use label encoding for ordinal data."

### Why Encode Categorical Variables?
- Many algorithms (K-means, SVM, neural networks) require numeric input
- Random Forest can handle factors, but some algorithms can't

### Label Encoding Example
```r
covid_imputed$cough_severity_encoded <- case_when(
  is.na(covid_imputed$cough_severity) | covid_imputed$cough_severity == "" ~ 0,
  covid_imputed$cough_severity == "Mild" ~ 1,
  covid_imputed$cough_severity == "Moderate" ~ 2,
  covid_imputed$cough_severity == "Severe" ~ 3,
  TRUE ~ 0
)
```

**Explain `case_when()` Function:**
- **Purpose:** Vectorized if-else statements (like SQL CASE WHEN)
- **Syntax:** `condition ~ value` pairs
- **Evaluation:** Checks conditions in order, returns first matching value
- **`TRUE ~ 0`:** Default case (catches anything not matched above)

**Why This Encoding?**
- **Ordinal data:** The categories have a natural order (None < Mild < Moderate < Severe)
- **Numeric mapping preserves order:** 0 < 1 < 2 < 3
- This is appropriate for ordinal data

**Alternative Encoding Methods to Mention:**
- **One-hot encoding:** For nominal (non-ordered) categories
  - Creates binary columns for each category
  - Example: "Red", "Blue", "Green" → 3 columns
- **Target encoding:** Uses target variable statistics (advanced)

**Verify the Encoding:**
```r
table(covid_imputed$cough_severity_encoded)
```
- Always check your transformations!
- Should see counts for 0, 1, 2, 3

---

## Section 5: Feature Scaling (10 minutes)

**What to Say:**
"Some algorithms use distance calculations. If features have very different scales, the larger-scaled features will dominate. We need to standardize features so they're on the same scale."

### Why Scaling Matters
**Demonstrate the Problem:**
```r
cat("Temperature range:", range(covid_imputed$temperature, na.rm = TRUE), "\n")
```

**Explain:**
- Temperature might range from 36-40°C
- Binary symptoms are 0 or 1
- In distance calculations (like K-means), temperature differences (e.g., 4 units) will dominate over symptom differences (0 or 1)
- This can make the algorithm ignore important features!

### Standardization (Z-score Normalization)
```r
temperature_scaled <- scale(covid_imputed$temperature)
```

**Explain `scale()` Function:**
- **Purpose:** Standardizes a vector to mean=0, standard deviation=1
- **Formula:** `(x - mean(x)) / sd(x)`
- **Parameters:**
  - `x`: The vector to scale
  - `center = TRUE` (default): Subtract the mean
  - `scale = TRUE` (default): Divide by standard deviation
- **Returns:** Scaled vector with mean ≈ 0, sd ≈ 1

**Verify the Scaling:**
```r
cat("Original - Mean:", round(mean(covid_imputed$temperature), 2), 
    "SD:", round(sd(covid_imputed$temperature), 2), "\n",
    "Scaled - Mean:", round(mean(temperature_scaled), 4), 
    "SD:", round(sd(temperature_scaled), 4), "\n")
```

**Key Points:**
- Original data: mean ≈ 37, SD ≈ 1
- Scaled data: mean ≈ 0, SD ≈ 1
- All features now on same scale!

**When to Scale:**
- ✅ **Required for:** K-means, SVM, KNN, neural networks, PCA
- ❌ **Not needed for:** Random Forest, decision trees (scale-invariant)
- ❌ **Not needed for:** Logistic regression (but can help convergence)

**Alternative Scaling Methods:**
- **Min-Max scaling:** Scales to [0, 1] range
- **Robust scaling:** Uses median and IQR (better for outliers)

**Important Note:**
- We demonstrate scaling here, but will apply it to all features in Part 2 when we use K-means
- Random Forest (Part 3) doesn't need scaling

---

## Section 6: Prepare Final Dataset (8 minutes)

**What to Say:**
"Now we'll create our final clean dataset with all the features we want to use for machine learning. We'll select only the relevant columns and ensure everything is ready."

### Select Features
```r
covid_ml <- covid_imputed %>%
  select(
    # Target
    covid19_test_results,
    # Numeric features
    temperature,
    # Binary features
    high_risk_exposure_occupation,
    high_risk_interactions,
    labored_respiration,
    # ... (other features)
    # Encoded categorical
    cough_severity_encoded
  )
```

**Explain the `select()` Function:**
- **Purpose:** Selects specific columns from a dataframe
- **Parameters:** Column names (unquoted) or column indices
- **Pipe operator `%>%`:** Passes the left-hand side as first argument to right-hand side
  - `x %>% f(y)` is equivalent to `f(x, y)`
  - Makes code more readable

**Why Select Specific Columns?**
- Remove irrelevant features
- Remove original categorical columns (we have encoded versions)
- Keep only what we need for ML

### Verify Final Dataset
```r
cat("Final dataset dimensions:", dim(covid_ml), "\n",
    "Missing values:", sum(is.na(covid_ml)), "\n")
```

**Explain `dim()` Function:**
- Returns: `[number of rows, number of columns]`
- Quick check of dataset size

**Final Checks:**
- ✅ No missing values
- ✅ All features are numeric (or factors for Random Forest)
- ✅ Target variable included

### Save Cleaned Data
```r
write.csv(covid_ml, "data/covid_ml_clean.csv", row.names = FALSE)
```

**Explain `write.csv()` Function:**
- **Parameters:**
  - `x`: The dataframe to write
  - `file`: Output file path
  - `row.names = FALSE`: Don't write row numbers (usually not needed)
- **Why save?** So we can use this cleaned data in Parts 2 and 3 without re-running preprocessing

---

## Section 7: Summary and Transition (5 minutes)

**What to Say:**
"Great job! We've completed the data preprocessing. Let's recap what we learned and prepare for the next parts."

### Key Takeaways
1. **Always explore your data first** - understand what you're working with
2. **Handle missing values** - remove or impute (we used median imputation)
3. **Encode categorical variables** - convert to numeric when needed
4. **Scale features** - when using distance-based algorithms

### What's Next?
- **Part 2:** K-means clustering (unsupervised learning) - will use scaling
- **Part 3:** Random Forest (supervised learning) - doesn't need scaling

**Questions to Ask:**
- "Why do you think we need to scale for K-means but not Random Forest?"
- "What would happen if we didn't handle missing values?"
- "When might you choose to remove rows instead of imputing?"

---

## Common Student Questions & Answers

**Q: Why median instead of mean for imputation?**
A: Median is robust to outliers. If you have a few extreme values, mean would be skewed, but median represents the "typical" value better.

**Q: Can we use the same scaling for training and test data?**
A: No! You must fit the scaling on training data only, then apply those same parameters to test data. Otherwise, you're "cheating" by using information from the test set.

**Q: What if we have too many missing values?**
A: If >50% of a column is missing, consider removing that feature entirely. If >50% of a row is missing, consider removing that observation.

**Q: Do we always need to encode categorical variables?**
A: It depends on the algorithm. Random Forest can handle factors directly, but K-means requires numeric input.

---

## Teaching Tips

1. **Emphasize the "why"** - Students often skip preprocessing, so explain why each step matters
2. **Show the data at each step** - Use `head()`, `str()`, or `View()` to show transformations
3. **Common mistakes to warn about:**
   - Forgetting to handle missing values
   - Scaling test data independently
   - Encoding categorical variables incorrectly (losing order information)
4. **Interactive elements:**
   - Ask students to check their data at each step
   - Have them identify which columns need what type of preprocessing
   - Encourage experimentation with different imputation methods
