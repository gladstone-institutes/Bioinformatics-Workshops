# Teaching Narrative: Part 3 - Random Forest

## Introduction (8 minutes)

**Opening Statement:**
"Welcome to Part 3! Now we're moving from unsupervised learning to supervised learning. We'll use Random Forest - a powerful and popular machine learning algorithm that can predict COVID test results from symptoms."

**Key Concepts to Introduce:**

1. **Supervised Learning:**
   - We have a target variable (what we want to predict)
   - We train on labeled data (known outcomes)
   - Goal: Learn patterns to predict new cases

2. **Random Forest:**
   - **Ensemble method:** Combines many decision trees
   - **Two sources of randomness:**
     - Bootstrap sampling (different data for each tree)
     - Feature sampling (different features considered at each split)
   - **Why it works:** Many diverse trees → better predictions than one tree

3. **Decision Trees (Brief Intro):**
   - Simple if-then rules
   - Example: "If temperature > 38°C AND cough = Yes, then likely Positive"
   - Easy to interpret, but can overfit

**Real-world Applications:**
- Medical diagnosis
- Credit risk assessment
- Image classification
- Feature importance analysis

**Why Random Forest?**
- ✅ Handles mixed data types (numeric + categorical)
- ✅ Provides feature importance
- ✅ Robust to overfitting
- ✅ No scaling required!
- ✅ Works well out-of-the-box

---

## Section 1: Load Libraries and Data (3 minutes)

**What to Say:**
"Let's load our libraries and the cleaned dataset from Part 1."

### Load Libraries

```r
library(tidyverse)    # Data manipulation and visualization
library(caret)        # Classification and regression training
library(randomForest) # Random Forest implementation
library(ggplot2)      # Plotting
```

**Explain Each Library:**

- **`tidyverse`:** Collection of packages for data science
  - `dplyr` for data manipulation
  - `ggplot2` for visualization

- **`caret`:** "Classification And REgression Training"
  - Provides unified interface for many ML algorithms
  - Handles cross-validation, hyperparameter tuning
  - Functions: `train()`, `trainControl()`, `createDataPartition()`

- **`randomForest`:** Classic Random Forest implementation
  - Fast, well-tested
  - Functions: `randomForest()`, `importance()`, `varImpPlot()`

- **`ggplot2`:** Grammar of graphics for creating plots

### Set Seed

```r
set.seed(123)
```

**Explain:**
- Random Forest uses randomness (bootstrap sampling, feature sampling)
- Setting seed ensures reproducible results
- **Critical for teaching:** Everyone gets same results!

### Load and Prepare Data

```r
covid_ml <- read.csv("data/covid_ml_clean.csv")
covid_ml$covid19_test_results <- as.factor(covid_ml$covid19_test_results)
```

**Explain `as.factor()`:**
- **Purpose:** Converts a vector to a factor (categorical variable)
- **Why needed:** Random Forest for classification requires target to be factor
- **What it does:** Treats values as categories, not numbers
- Levels: "Negative", "Positive" (or whatever your classes are)

**Teaching Point:**
"Unlike K-means where we ignored the target, here the target is central - it's what we're trying to predict!"

---

## Section 2: Train-Test Split (10 minutes)

**What to Say:**
"Before training any model, we need to split our data. This is crucial for evaluating how well our model will work on new, unseen data."

### Why Split the Data?

**The Problem:**
- If we train and test on the same data, the model might "memorize" the training data
- This gives overly optimistic performance estimates
- Model might not generalize to new data

**The Solution:**
- **Training set (70%):** Used to train/fit the model
- **Test set (30%):** Used to evaluate performance on unseen data
- Test set is "held out" - never used during training

**Analogy:**
"Think of it like studying for an exam. If you only test yourself on the exact problems you studied, you'll do great - but that doesn't mean you understand the material. The test set is like a real exam with new problems."

### Perform the Split

```r
train_index <- createDataPartition(
  covid_ml$covid19_test_results, 
  p = 0.7, 
  list = FALSE
)
```

**Explain `createDataPartition()` Function:**

**Purpose:** Creates indices for splitting data into training and testing sets

**Parameters:**

1. **`y` (first argument):**
   - The target variable (vector or factor)
   - Used to ensure balanced split across classes

2. **`p = 0.7`:**
   - Proportion of data for training set
   - 0.7 = 70% training, 30% testing
   - Common splits: 70/30, 80/20, or 90/10 (for large datasets)

3. **`list = FALSE`:**
   - Returns a vector of indices (TRUE) or matrix (FALSE)
   - `FALSE` gives a vector, easier to use

**Why Stratified Split?**
- `createDataPartition()` maintains class proportions
- If 50% Positive in full data → ~50% Positive in training AND test sets
- Prevents class imbalance issues

**Returns:** Vector of row indices for training set

### Create Training and Test Sets

```r
train_data <- covid_ml[train_index, ]
test_data <- covid_ml[-train_index, ]
```

**Explain Indexing:**

- **`covid_ml[train_index, ]`:** 
  - Selects rows in `train_index`
  - `, ]` means "all columns"
  - Creates training set

- **`covid_ml[-train_index, ]`:**
  - `-` negates the indices
  - Selects rows NOT in `train_index`
  - Creates test set

**Verify the Split:**

```r
cat("Training set size:", nrow(train_data), "\n",
    "Test set size:", nrow(test_data), "\n")
```

**Teaching Point:**
"Always verify your split! Check sizes, check class distributions in both sets."

### Important Note: No Scaling Needed!

**What to Say:**
"Unlike K-means, Random Forest doesn't require feature scaling. Tree-based methods are invariant to monotonic transformations - they only care about the order of values, not their scale."

**Why?**
- Decision trees split on thresholds (e.g., "temperature > 38")
- The threshold value matters, but the scale doesn't
- Whether temperature is 36-40 or 0.36-0.40, the splits work the same

**Teaching Point:**
"This is a major advantage of tree-based methods - less preprocessing needed!"

---

## Section 3: Train a Simple Random Forest (12 minutes)

**What to Say:**
"Let's start with a simple Random Forest to understand the basics. We'll use default parameters first, then tune them later."

### Train the Model

```r
rf_simple <- randomForest(
  covid19_test_results ~ .,  # Formula: predict target using all features
  data = train_data,
  ntree = 100,               # Number of trees
  importance = TRUE          # Calculate variable importance
)
```

**Explain `randomForest()` Function in Detail:**

**Purpose:** Trains a Random Forest model for classification or regression

**Parameters:**

1. **Formula Syntax: `y ~ .`**
   - **`y`:** Target variable (what we predict)
   - **`~`:** "is modeled by" or "depends on"
   - **`.`:** All other columns in the data
   - **Alternative:** `y ~ x1 + x2 + x3` (specify features explicitly)

2. **`data`:**
   - Training dataframe
   - Must contain target and all features in formula

3. **`ntree = 100`:**
   - Number of trees to grow
   - More trees = better predictions (but diminishing returns)
   - Default: 500 (we use 100 for speed)
   - **Trade-off:** More trees = slower training, but more stable predictions

4. **`importance = TRUE`:**
   - Calculate variable importance metrics
   - Needed for `importance()` and `varImpPlot()` functions
   - Slightly slower, but very useful!

**Other Important Parameters (Not Shown):**

- **`mtry`:** Number of features to consider at each split
  - Default: `sqrt(p)` for classification (p = number of features)
  - Controls randomness - lower mtry = more diverse trees
  - We'll tune this later!

- **`nodesize`:** Minimum size of terminal nodes
  - Default: 1 for classification
  - Larger = simpler trees, less overfitting

- **`maxnodes`:** Maximum number of terminal nodes
  - Default: unlimited
  - Controls tree complexity

**Returns:** Random Forest object containing:
- `$predicted`: Predictions on training data
- `$importance`: Variable importance matrix
- `$confusion`: Confusion matrix on training data
- `$err.rate`: Error rate as trees are added
- `$votes`: Class votes for each observation
- And more...

### View Model Summary

```r
print(rf_simple)
```

**What the Output Shows:**
- Type of random forest (classification)
- Number of trees
- No. of variables tried at each split (mtry)
- OOB estimate of error rate
- Confusion matrix

### Out-of-Bag (OOB) Error

**What to Say:**
"Random Forest has a built-in validation mechanism called Out-of-Bag (OOB) error. This is one of its best features!"

**How It Works:**

1. **Bootstrap Sampling:**
   - Each tree uses ~63% of data (sampled with replacement)
   - Remaining ~37% is "out-of-bag" (OOB) for that tree

2. **OOB Predictions:**
   - For each observation, predict using trees where it was OOB
   - Average these predictions

3. **OOB Error:**
   - Compare OOB predictions to true labels
   - Gives unbiased estimate of test error
   - **No separate validation set needed!**

**Why This Matters:**
- Provides performance estimate without using test set
- Uses all data for training (each observation is in some trees)
- More efficient than traditional train/validation split

**Teaching Point:**
"This is like having a built-in cross-validation! We can estimate performance without touching our test set."

---

## Section 4: Tune Random Forest with Cross-Validation (20 minutes)

**What to Say:**
"Now let's tune our hyperparameters to improve performance. We'll use cross-validation to find the best settings."

### Key Hyperparameters

**What to Say:**
"Random Forest has several hyperparameters, but the most important to tune is `mtry` - the number of features considered at each split."

**Hyperparameters Explained:**

1. **`ntree` (Number of Trees):**
   - More trees = better (but diminishing returns)
   - Usually 100-500 is sufficient
   - Easy to set, less critical to tune

2. **`mtry` (Features per Split):**
   - **Most important to tune!**
   - Default: `sqrt(p)` for classification
   - Lower mtry = more randomness = more diverse trees
   - Higher mtry = trees more similar to each other
   - **Sweet spot:** Usually between `sqrt(p)` and `p/3`

3. **`nodesize` and `maxnodes`:**
   - Control tree complexity
   - Usually left at defaults for first pass

### Set Up Cross-Validation

```r
cv_control <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)
```

**Explain `trainControl()` Function:**

**Purpose:** Defines how model training and resampling will be done

**Parameters:**

1. **`method = "cv"`:**
   - Cross-validation method
   - Options: `"cv"` (k-fold), `"repeatedcv"` (repeated k-fold), `"boot"` (bootstrap), etc.

2. **`number = 10`:**
   - Number of folds (for k-fold CV)
   - 10-fold CV: Split data into 10 parts
   - Train on 9 parts, validate on 1 part
   - Repeat 10 times, average results
   - **Why 10?** Good balance of bias and variance

3. **`classProbs = TRUE`:**
   - Calculate class probabilities (not just predictions)
   - Needed for ROC curve, probability-based metrics

4. **`summaryFunction = twoClassSummary`:**
   - Function to compute performance metrics
   - `twoClassSummary` calculates: Sensitivity, Specificity, ROC-AUC
   - For multi-class: use `multiClassSummary`

**How Cross-Validation Works:**

1. Split training data into 10 folds
2. For each fold:
   - Train on 9 folds
   - Validate on 1 fold
   - Calculate performance metrics
3. Average metrics across all 10 folds
4. This gives robust estimate of model performance

**Teaching Point:**
"Cross-validation uses only the training set - our test set is still untouched! This is proper ML practice."

### Define Tuning Grid

```r
mtry_grid <- expand.grid(mtry = c(2, 4, 6, 8, 10))
```

**Explain `expand.grid()`:**
- Creates a dataframe of all combinations of parameters
- Here: just mtry values to try
- For multiple parameters: creates full grid

**Why These Values?**
- `sqrt(p)` ≈ 4-5 (typical default)
- Testing range around default: 2, 4, 6, 8, 10
- Can expand range if needed

### Train with Cross-Validation

```r
rf_cv <- train(
  covid19_test_results ~ .,
  data = train_data,
  method = "rf",
  trControl = cv_control,
  tuneGrid = mtry_grid,
  ntree = 200,
  metric = "ROC"
)
```

**Explain `train()` Function:**

**Purpose:** Trains a model with cross-validation and hyperparameter tuning

**Parameters:**

1. **Formula:** Same as `randomForest()` - `y ~ .`

2. **`data`:** Training dataframe

3. **`method = "rf"`:**
   - Algorithm to use
   - `"rf"` = Random Forest
   - Caret supports 200+ algorithms!

4. **`trControl = cv_control`:**
   - Resampling method (from `trainControl()`)

5. **`tuneGrid = mtry_grid`:**
   - Grid of hyperparameters to try
   - Caret will test each combination
   - Returns best based on metric

6. **`ntree = 200`:**
   - Number of trees (passed to `randomForest()`)
   - More than simple model for better performance

7. **`metric = "ROC"`:**
   - Metric to optimize
   - `"ROC"` = Area Under ROC Curve (AUC-ROC)
   - Alternatives: `"Accuracy"`, `"Kappa"`

**What Happens:**
- For each mtry value (2, 4, 6, 8, 10):
  - Performs 10-fold CV
  - Calculates average ROC-AUC
- Selects mtry with highest ROC-AUC
- Trains final model on full training set with best mtry

### View Results

```r
rf_cv
```

**Output Shows:**
- Best mtry value
- Performance metrics for each mtry tried
- Final model details

### Plot Tuning Results

```r
plot(rf_cv, main = "Random Forest: Tuning mtry")
```

**What the Plot Shows:**
- X-axis: mtry values
- Y-axis: Performance metric (ROC-AUC)
- Helps visualize which mtry is best
- Shows if performance plateaus

**Teaching Point:**
"If the curve is flat, multiple mtry values work similarly. If there's a clear peak, that's the optimal value."

### Get Best mtry

```r
cat("Best mtry:", rf_cv$bestTune$mtry, "\n")
```

**Explain:**
- `rf_cv$bestTune` contains best hyperparameters
- `$mtry` gives the optimal mtry value
- Use this for final model or it's already used!

---

## Section 5: Variable Importance (15 minutes)

**What to Say:**
"One of Random Forest's best features is built-in variable importance. This tells us which features are most important for predictions - very useful for understanding your data!"

### Get Importance Metrics

```r
importance_df <- as.data.frame(importance(rf_simple))
importance_df$Variable <- rownames(importance_df)
head(importance_df)
```

**Explain `importance()` Function:**

**Purpose:** Extracts variable importance from Random Forest model

**Parameters:**
- `x`: Random Forest object
- `type`: Type of importance (1 = mean decrease accuracy, 2 = mean decrease Gini)
- `scale`: Whether to scale importance values

**Returns:** Matrix with importance metrics:
- **MeanDecreaseAccuracy:** How much accuracy decreases when variable is permuted
- **MeanDecreaseGini:** How much node purity (Gini) decreases without variable

**Mean Decrease Gini Explained:**
- **Gini impurity:** Measures how "pure" a node is
  - Pure node: all observations same class (Gini = 0)
  - Impure node: mixed classes (Gini > 0)
- **How it's calculated:**
  - For each tree, sum Gini decreases from splits using this variable
  - Average across all trees
- **Interpretation:** Higher = more important for classification

**Why Two Metrics?**
- **MeanDecreaseGini:** Faster to calculate, default
- **MeanDecreaseAccuracy:** More interpretable (directly measures prediction impact)

### Plot Variable Importance

```r
varImpPlot(
  rf_simple, 
  type = 2,  # Mean Decrease Gini
  main = "Variable Importance - Random Forest",
  n.var = 15
)
```

**Explain `varImpPlot()` Function:**

**Parameters:**

1. **`x`:** Random Forest object

2. **`type = 2`:**
   - `1` = Mean Decrease Accuracy
   - `2` = Mean Decrease Gini (default)

3. **`main`:** Plot title

4. **`n.var = 15`:**
   - Number of top variables to show
   - Shows most important features

**What the Plot Shows:**
- Bar plot of variable importance
- Sorted by importance (highest at top)
- Helps identify which symptoms/features matter most

**Teaching Points:**
- "This is actionable insight - which symptoms should doctors focus on?"
- "High importance doesn't mean causation - just strong association"
- "Use this to simplify models (remove low-importance features)"

---

## Section 6: Make Predictions (12 minutes)

**What to Say:**
"Now let's use our trained model to make predictions on the test set. This is where we see how well our model generalizes to new data."

### Predict Classes

```r
rf_predictions <- predict(rf_cv, newdata = test_data)
head(rf_predictions)
```

**Explain `predict()` Function:**

**Purpose:** Makes predictions on new data

**Parameters:**

1. **`object`:** Trained model (from `train()` or `randomForest()`)

2. **`newdata`:** New dataframe with same features as training data
   - Must have same column names and types
   - Missing features will cause errors

3. **`type`:** Type of prediction
   - Default: class predictions (for classification)
   - `type = "prob"`: class probabilities
   - `type = "vote"`: raw votes from trees

**Returns:** Vector of predicted classes (same length as test set)

**Teaching Point:**
"These are predictions on data the model has never seen - this is the real test!"

### Predict Probabilities

```r
rf_proba <- predict(rf_cv, newdata = test_data, type = "prob")
head(rf_proba)
```

**Explain:**
- `type = "prob"` returns probabilities for each class
- Returns dataframe with columns for each class
- Each row sums to 1.0
- Probabilities show model's confidence

**Why Probabilities Matter:**
- More informative than just class predictions
- Can adjust decision threshold
- Useful for risk stratification

### Custom Threshold

```r
threshold <- 0.3
rf_predictions_custom <- ifelse(
  rf_proba$Positive > threshold, 
  "Positive", 
  "Negative"
)
rf_predictions_custom <- factor(
  rf_predictions_custom, 
  levels = c("Negative", "Positive")
)
```

**Explain the Concept:**

**Default Threshold (0.5):**
- If P(Positive) > 0.5 → predict Positive
- Balanced approach

**Custom Threshold:**
- **Lower threshold (e.g., 0.3):**
  - More sensitive (catches more positives)
  - But more false positives
  - Use when missing a positive case is costly (e.g., disease screening)

- **Higher threshold (e.g., 0.7):**
  - More specific (fewer false positives)
  - But may miss some positives
  - Use when false alarms are costly

**The Trade-off:**
- **Sensitivity vs Specificity**
- Can't maximize both simultaneously
- Choose based on your application's needs

**Explain the Code:**
- `ifelse(condition, value_if_true, value_if_false)`: Vectorized if-else
- `factor(..., levels = ...)`: Ensures factor has correct levels
- Compare counts to see impact of threshold change

**Teaching Point:**
"In medical diagnosis, you might prefer high sensitivity (catch all cases) even with more false alarms. In spam detection, you might prefer high specificity (fewer false alarms)."

---

## Section 7: Evaluate Model Performance (20 minutes)

**What to Say:**
"Now for the moment of truth - how well does our model actually perform? We'll use a confusion matrix and various metrics to evaluate."

### Confusion Matrix Concept

**What to Say:**
"A confusion matrix is a table that shows how well our model classified each class. It's the foundation for all classification metrics."

**Structure:**
```
                    Predicted
                 Negative  Positive
Actual Negative    TN        FP
       Positive    FN        TP
```

**Explain Each Cell:**

- **TN (True Negative):** Correctly predicted negative
- **TP (True Positive):** Correctly predicted positive
- **FN (False Negative):** Missed a positive case (Type II error)
- **FP (False Positive):** False alarm (Type I error)

**Teaching Point:**
"Think about the cost of each error type. In medical diagnosis, false negatives (missing a disease) are often worse than false positives (unnecessary test)."

### Create Confusion Matrix

```r
cm_rf <- confusionMatrix(
  rf_predictions, 
  test_data$covid19_test_results, 
  positive = "Positive"
)
print(cm_rf)
```

**Explain `confusionMatrix()` Function:**

**Purpose:** Creates confusion matrix and calculates performance metrics

**Parameters:**

1. **`data`:** Vector of predictions

2. **`reference`:** Vector of true labels
   - Must have same length as predictions

3. **`positive`:** Which class is "positive"
   - Needed to calculate sensitivity, specificity correctly
   - Here: "Positive" is the positive class

**Returns:** Confusion matrix object with:
- `$table`: The confusion matrix table
- `$overall`: Overall metrics (Accuracy, Kappa, etc.)
- `$byClass`: Class-specific metrics (Sensitivity, Specificity, etc.)

### Understanding the Metrics

```r
accuracy <- cm_rf$overall["Accuracy"]
sensitivity <- cm_rf$byClass["Sensitivity"]
specificity <- cm_rf$byClass["Specificity"]
```

**Explain Each Metric:**

1. **Accuracy:**
   - Formula: `(TP + TN) / (TP + TN + FP + FN)`
   - Proportion of correct predictions
   - **Limitation:** Can be misleading with imbalanced classes
   - Example: 95% accuracy sounds great, but if 95% are negative, always predicting negative gives 95% accuracy!

2. **Sensitivity (Recall, True Positive Rate):**
   - Formula: `TP / (TP + FN)`
   - Of all actual positives, how many did we catch?
   - **Interpretation:** Ability to detect positive cases
   - **When important:** Medical screening, fraud detection
   - **Also called:** Recall, Hit Rate

3. **Specificity (True Negative Rate):**
   - Formula: `TN / (TN + FP)`
   - Of all actual negatives, how many did we correctly identify?
   - **Interpretation:** Ability to avoid false alarms
   - **When important:** When false positives are costly

**Other Important Metrics:**

- **Precision (Positive Predictive Value):**
  - `TP / (TP + FP)`
  - Of all predicted positives, how many were actually positive?
  - Important when false positives are costly

- **F1-Score:**
  - Harmonic mean of Precision and Sensitivity
  - `2 * (Precision * Sensitivity) / (Precision + Sensitivity)`
  - Balances precision and recall

- **ROC-AUC:**
  - Area Under ROC Curve
  - Measures ability to distinguish between classes
  - Range: 0.5 (random) to 1.0 (perfect)
  - Threshold-independent metric

**Teaching Points:**
- "No single metric tells the whole story - look at multiple metrics"
- "Choose metrics based on your application's needs"
- "Accuracy can be misleading with imbalanced data"

### Visualize Confusion Matrix

```r
cm_df <- as.data.frame(cm_rf$table)

ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), size = 8, color = "white") +
  scale_fill_gradient(low = "#27ae60", high = "#c0392b") +
  labs(title = "Confusion Matrix - Random Forest",
       x = "Actual",
       y = "Predicted") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12))
```

**Explain the Visualization:**

- **Heatmap:** Color intensity shows count
- **Diagonal (green):** Correct predictions (TN, TP)
- **Off-diagonal (red):** Errors (FP, FN)
- **Text labels:** Actual counts in each cell

**How to Read:**
- Dark green diagonal = good performance
- Red off-diagonal = errors to investigate
- Compare sizes to see which errors are more common

**Teaching Point:**
"Visualization helps quickly see where the model struggles. Are there more false positives or false negatives?"

---

## Section 8: Summary and Wrap-up (10 minutes)

**What to Say:**
"Congratulations! You've completed all three parts of the workshop. Let's recap what we learned about Random Forest and machine learning in general."

### What We Learned About Random Forest

1. **Ensemble Method:**
   - Combines many decision trees
   - Bootstrap + feature sampling create diversity
   - Reduces overfitting compared to single tree

2. **No Scaling Required:**
   - Tree-based methods are scale-invariant
   - Major advantage over K-means, SVM, neural networks

3. **Built-in Features:**
   - Variable importance (understand which features matter)
   - OOB error (validation without separate set)
   - Handles mixed data types

4. **Hyperparameter Tuning:**
   - Most important: `mtry` (features per split)
   - Use cross-validation to find optimal values
   - More trees usually better (but diminishing returns)

5. **Evaluation:**
   - Confusion matrix shows detailed performance
   - Multiple metrics (accuracy, sensitivity, specificity)
   - Choose metrics based on application needs

### When to Use Random Forest

**✅ Good For:**
- Mixed numeric and categorical features
- Want to understand feature importance
- Large datasets
- Don't want to tune many hyperparameters
- Need robust baseline model
- Non-linear relationships

**❌ Not Ideal For:**
- Very high-dimensional data (p >> n)
- Need interpretable individual predictions
- Extremely large datasets (might be slow)
- Linear relationships (linear models might be simpler)

### Comparison: K-means vs Random Forest

| Aspect | K-means | Random Forest |
|--------|---------|---------------|
| Type | Unsupervised | Supervised |
| Target Variable | Not used | Required |
| Scaling | Required | Not needed |
| Output | Clusters | Predictions + Importance |
| Interpretability | Moderate | High (importance) |
| Use Case | Finding patterns | Making predictions |

### Workshop Complete!

**What You've Learned:**
1. ✅ **Data Preprocessing:** Handle missing values, encode categories, scale features
2. ✅ **Unsupervised Learning:** K-means clustering to find patterns
3. ✅ **Supervised Learning:** Random Forest for classification
4. ✅ **Model Evaluation:** Confusion matrices, metrics, visualization

**Next Steps:**
- Try other algorithms (SVM, neural networks, XGBoost)
- Explore other evaluation metrics (ROC curves, precision-recall curves)
- Learn about feature engineering
- Practice on your own datasets!

---

## Common Student Questions & Answers

**Q: Why do we need both training and test sets?**
A: The training set teaches the model, the test set evaluates it. If we evaluate on training data, we get overly optimistic results (the model "memorized" the data). The test set tells us how well the model will work on truly new data.

**Q: What's the difference between OOB error and test error?**
A: OOB error is calculated during training using bootstrap samples. Test error uses completely held-out data. Both estimate generalization error, but test error is the "gold standard" final evaluation.

**Q: How many trees should I use?**
A: Usually 100-500 is sufficient. More trees = better (but diminishing returns) and slower. Start with 100-200, increase if needed. The OOB error stabilizes as you add trees.

**Q: What if my classes are imbalanced?**
A: Random Forest can handle imbalance, but you might want to:
- Use `classwt` parameter to weight classes
- Use different evaluation metrics (F1, ROC-AUC instead of accuracy)
- Consider resampling techniques (SMOTE, etc.)

**Q: Can Random Forest overfit?**
A: Yes, but less than single decision trees. Signs of overfitting:
- Training accuracy much higher than test accuracy
- Very deep trees (high maxnodes)
- Too many trees with very low mtry

**Q: How do I interpret variable importance?**
A: Higher values = more important for predictions. But:
- Importance doesn't mean causation
- Correlated features might have lower individual importance
- Use domain knowledge to validate importance rankings

**Q: What's the difference between mtry and ntree?**
A: 
- **`ntree`:** Number of trees in the forest (more = better, but slower)
- **`mtry`:** Number of features considered at each split (controls diversity of trees)

**Q: When should I use probabilities vs class predictions?**
A: Use probabilities when:
- You need to adjust decision threshold
- You want to rank observations by risk
- You're building a larger system that uses confidence scores

Use class predictions when:
- You just need yes/no answers
- Default threshold (0.5) is appropriate

---

## Teaching Tips

1. **Emphasize the supervised nature:**
   - Constantly contrast with K-means (which ignored target)
   - Show how target variable guides learning

2. **Make evaluation concrete:**
   - Use confusion matrix to show exactly what "good" and "bad" predictions look like
   - Connect metrics to real-world consequences

3. **Highlight Random Forest advantages:**
   - No scaling needed (contrast with K-means)
   - Built-in importance (very useful!)
   - Robust to overfitting

4. **Common mistakes to warn about:**
   - Evaluating on training data (show the difference!)
   - Not setting seed (results won't be reproducible)
   - Ignoring class imbalance (accuracy can be misleading)

5. **Interactive elements:**
   - Have students try different mtry values and compare
   - Ask them to interpret variable importance in biological context
   - Discuss what threshold they'd use for medical diagnosis

6. **Connect to biology:**
   - "Which symptoms are most predictive? Does this match clinical knowledge?"
   - "How would you use this model in practice?"
   - "What are the ethical considerations of ML in healthcare?"

7. **Build intuition:**
   - Explain decision trees with simple examples
   - Show how ensemble reduces variance
   - Use analogies (committee voting, etc.)
