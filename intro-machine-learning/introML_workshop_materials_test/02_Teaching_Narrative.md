# Teaching Narrative: Part 2 - K-means Clustering

## Introduction (5 minutes)

**Opening Statement:**
"Welcome to Part 2! Now that we have clean, preprocessed data, we're going to learn about K-means clustering - an unsupervised learning method. Unlike supervised learning where we predict a known outcome, clustering finds hidden patterns in data without knowing the 'answer' beforehand."

**Key Concepts to Introduce:**
- **Unsupervised learning:** No target variable, finding patterns in data
- **Clustering:** Grouping similar observations together
- **K-means:** A specific clustering algorithm that partitions data into K groups

**Real-world Applications:**
- Customer segmentation (marketing)
- Image compression
- Gene expression analysis
- Anomaly detection

---

## Section 1: Load Libraries and Data (3 minutes)

**What to Say:**
"Let's load the libraries we need and import our cleaned dataset from Part 1."

### Load Libraries
```r
library(cluster)      # Clustering algorithms and utilities
library(factoextra)   # Visualization tools for clustering
library(ggplot2)      # Plotting
library(dplyr)        # Data manipulation
library(SamSPECTRAL)  # For automated elbow detection
```

**Explain Each Library:**
- **`cluster`:** Contains `kmeans()` function and other clustering methods
- **`factoextra`:** Provides `fviz_cluster()` and `fviz_nbclust()` for beautiful visualizations
- **`ggplot2`:** Grammar of graphics for creating plots
- **`dplyr`:** Data manipulation (select, filter, mutate, etc.)
- **`SamSPECTRAL`:** Contains `kneepointDetection()` for finding optimal K

### Set Seed
```r
set.seed(123)
```

**Explain `set.seed()`:**
- **Purpose:** Makes random processes reproducible
- **Parameter:** Any integer (the "seed" value)
- **Why important:** K-means uses random initialization, so results vary without a seed
- **Teaching point:** Always set seed when using random algorithms for reproducibility!

### Load Data
```r
covid_ml <- read.csv("data/covid_ml_clean.csv")
head(covid_ml)
```

**Explain:**
- Loading the cleaned dataset from Part 1
- `head()` shows first 6 rows - quick data check

---

## Section 2: Prepare Data for K-means (10 minutes)

**What to Say:**
"K-means has specific requirements. Let's make sure our data meets them."

### K-means Requirements

**⚠️ Three Critical Requirements:**

1. **No missing values** ✅ (We handled this in Part 1!)
2. **Numeric data only** ✅ (We encoded categorical variables!)
3. **Scaled features** ⚠️ (We need to do this now!)

**Why These Requirements?**
- K-means uses **Euclidean distance** to measure similarity
- Missing values break distance calculations
- Categorical variables can't be used in distance formulas
- Unscaled features cause large-scale features to dominate

### Step 1: Select Features (Exclude Target)

```r
cluster_data <- covid_ml %>%
  select(-covid19_test_results) %>% # Remove target
  as.data.frame()
```

**Explain:**
- **Why exclude target?** K-means is unsupervised - we don't use the target variable
- **`select(-column_name)`:** The minus sign removes that column
- **`as.data.frame()`:** Ensures we have a dataframe (some functions require this)

**Teaching Point:**
"Even though we have a target variable, we're ignoring it for clustering. This is the key difference between supervised and unsupervised learning!"

### Step 2: Scale the Data

```r
cluster_scaled <- scale(cluster_data)
```

**Explain `scale()` Function in Detail:**

**What it does:**
- Standardizes each column to mean=0, standard deviation=1
- Formula: `(x - mean(x)) / sd(x)`

**Parameters:**
- `x`: Matrix or dataframe to scale
- `center = TRUE` (default): Subtract column means
- `scale = TRUE` (default): Divide by column standard deviations

**Returns:**
- Matrix with same dimensions
- Each column has mean ≈ 0, sd ≈ 1

**Why scale ALL features together?**
- K-means calculates distances across all features
- If temperature (36-40) isn't scaled but symptoms (0-1) are, temperature dominates
- Scaling puts everything on equal footing

### Step 3: Verify Scaling

```r
cat("Scaled data - Column means (should be ~0):\n")
round(colMeans(cluster_scaled), 10)[1:5]

cat("\nScaled data - Column SDs (should be ~1):\n")
round(apply(cluster_scaled, 2, sd), 10)[1:5]
```

**Explain the Verification:**

**`colMeans()` function:**
- Calculates mean of each column
- Returns a vector of means
- Should be approximately 0 for scaled data

**`apply()` function:**
- **Parameters:**
  - `X`: Matrix or dataframe
  - `MARGIN = 2`: Apply function to columns (1 = rows)
  - `FUN = sd`: Function to apply (standard deviation)
- **Returns:** Vector of standard deviations for each column
- Should be approximately 1 for scaled data

**`round(..., 10)`:**
- Rounds to 10 decimal places
- Shows very small numbers (essentially 0)

**Why verify?**
- Always check your preprocessing worked correctly!
- If means aren't ~0 or SDs aren't ~1, something went wrong

---

## Section 3: Choosing K - The Elbow Method (15 minutes)

**What to Say:**
"The hardest part of K-means is choosing K - how many clusters? There's no perfect answer, but we have methods to help us decide."

### The Challenge of Choosing K

**Key Points:**
- K must be specified before running the algorithm
- Too few clusters: groups are too broad, miss important patterns
- Too many clusters: overfitting, each point might be its own cluster
- There's no "correct" K - it depends on your goals

### The Elbow Method

**Concept:**
- Plot **Within-Cluster Sum of Squares (WSS)** vs K
- WSS measures how compact/tight clusters are
- Lower WSS = tighter clusters (better)
- As K increases, WSS always decreases
- Look for the "elbow" - where the rate of decrease slows

**Analogy:**
"Imagine you're wrapping rubber bands around groups of points. With K=1, you need one big band (high WSS). As K increases, bands get tighter (WSS decreases). The elbow is where adding more bands doesn't help much."

### Calculate WSS for Different K Values

```r
fviz.p <- fviz_nbclust(
  x = cluster_scaled, 
  FUNcluster = kmeans, 
  method = "wss", 
  k.max = 10
) +
  labs(
    title = "Elbow Method for Optimal K",
    subtitle = "Look for the 'elbow' where WSS decrease slows"
  ) +
  theme_minimal()
```

**Explain `fviz_nbclust()` Function:**

**Purpose:** Visualizes optimal number of clusters using various methods

**Parameters:**
- **`x`:** Scaled data matrix/dataframe
- **`FUNcluster`:** Clustering function to use (`kmeans`, `hclust`, etc.)
- **`method`:** Method for determining optimal K
  - `"wss"`: Within-cluster sum of squares (elbow method)
  - `"silhouette"`: Average silhouette width
  - `"gap_stat"`: Gap statistic
- **`k.max`:** Maximum K to test (tests K=1 to k.max)

**Returns:** A ggplot object that can be customized

**What the Plot Shows:**
- X-axis: Number of clusters (K)
- Y-axis: WSS (within-cluster sum of squares)
- Look for the "elbow" - sharp bend in the curve
- The elbow suggests optimal K (where adding more clusters doesn't help much)

**Teaching Tips:**
- Point out that the curve always decreases (more clusters = lower WSS)
- Emphasize looking for where the decrease rate changes
- Sometimes the elbow isn't obvious - that's okay!

### Automated Elbow Detection

```r
fviz.y <- ggplot_build(fviz.p)$data[[1]][,"y"]
detected <- kneepointDetection(vect = fviz.y)
cat("Detected optimal K (knee point):", detected$MinIndex, "\n")
```

**Explain This Code:**

**`ggplot_build()`:**
- Extracts underlying data from a ggplot object
- `$data[[1]]` gets the first layer's data
- `[,"y"]` extracts the y-values (WSS values)

**`kneepointDetection()` Function:**
- **Purpose:** Automatically finds the "knee" or "elbow" point
- **Parameter:** `vect` - vector of y-values (WSS values)
- **Returns:** List with `MinIndex` - the K value at the elbow

**Why Use This?**
- Provides objective suggestion for optimal K
- But remember: optimal K depends on your goals, not just statistics!

**Teaching Point:**
"Automated methods are helpful, but domain knowledge matters too. If you're looking for 3 customer segments, use K=3 even if the elbow suggests K=5!"

---

## Section 4: Train K-means (10 minutes)

**What to Say:**
"Now that we've chosen K (or at least have a suggestion), let's train our K-means model. We'll use K=2 as an example."

### Train the Model

```r
fit <- kmeans(
  cluster_scaled,
  centers = 2,
  nstart = 50,
  iter.max = 100,
  algorithm = "Lloyd"
)
```

**Explain `kmeans()` Function in Detail:**

**Purpose:** Performs K-means clustering on a dataset

**Parameters:**

1. **`x` (first argument):**
   - Scaled numeric data (matrix or dataframe)
   - Must be numeric, no missing values

2. **`centers`:**
   - Number of clusters (K) OR
   - Matrix of initial cluster centers
   - We use K=2 (integer)

3. **`nstart = 50`:**
   - Number of random initializations
   - K-means is sensitive to initial cluster centers (randomly chosen)
   - Algorithm runs 50 times with different starting points
   - Returns the best result (lowest WSS)
   - **Why important:** Prevents getting stuck in local optima
   - **Trade-off:** Higher nstart = better results but slower

4. **`iter.max = 100`:**
   - Maximum number of iterations per run
   - Algorithm stops if it converges before this
   - **What happens:** Algorithm iteratively updates cluster centers until convergence

5. **`algorithm = "Lloyd"`:**
   - Classic K-means algorithm (default)
   - Alternative: `"MacQueen"` (slightly different update rule)
   - Lloyd is most common

**Returns:** A list containing:
- `$cluster`: Vector of cluster assignments (1, 2, ..., K)
- `$centers`: Matrix of cluster centers (centroids)
- `$totss`: Total sum of squares
- `$withinss`: Within-cluster sum of squares for each cluster
- `$tot.withinss`: Total within-cluster sum of squares
- `$betweenss`: Between-cluster sum of squares
- `$size`: Number of points in each cluster
- `$iter`: Number of iterations until convergence

### View Cluster Assignments

```r
head(fit$cluster)
```

**Explain:**
- `fit$cluster` is a vector of cluster assignments
- Each observation gets assigned to cluster 1 or 2 (since K=2)
- Same length as number of rows in data
- Use `table(fit$cluster)` to see distribution

**Teaching Point:**
"Each patient is now assigned to a cluster. We don't know what these clusters represent yet - that's what makes it unsupervised!"

---

## Section 5: Visualize Clusters (8 minutes)

**What to Say:**
"Our data has many dimensions (features), but we can only visualize 2-3 dimensions. We'll use PCA (Principal Component Analysis) to reduce dimensions for visualization."

### Create Visualization

```r
fviz_cluster(
  fit, 
  data = cluster_scaled, 
  palette = c("#e74c3c", "#3498db"), 
  ggtheme = theme_minimal()
)
```

**Explain `fviz_cluster()` Function:**

**Purpose:** Creates a 2D visualization of clusters using PCA

**Parameters:**

1. **`object` (first argument):**
   - K-means result object (from `kmeans()`)

2. **`data`:**
   - Original scaled data used for clustering
   - Needed to project points onto PCA space

3. **`palette`:**
   - Colors for each cluster
   - Vector of color codes (hex, names, etc.)
   - Length should match number of clusters

4. **`ggtheme`:**
   - ggplot2 theme for styling
   - `theme_minimal()` gives clean, simple look

**What the Plot Shows:**
- Points colored by cluster assignment
- Ellipses around each cluster (confidence regions)
- Cluster centers (centroids) marked
- **Important:** This is a 2D projection - real clusters exist in high-dimensional space!

**How It Works:**
- Uses PCA to reduce dimensions to 2D
- Projects all points onto first 2 principal components
- Shows clusters in this reduced space
- **Limitation:** Some cluster separation might be lost in projection

**Teaching Points:**
- "The visualization is helpful but not perfect - it's a 2D view of high-dimensional data"
- "Points close together in the plot are similar across all features"
- "The ellipses show the spread of each cluster"

---

## Section 6: Compare Clusters with COVID Status (10 minutes)

**What to Say:**
"Even though K-means is unsupervised, we can check if the clusters we found correspond to something meaningful - like COVID test results. This is validation that our clusters make biological sense."

### Add Cluster Assignments to Data

```r
covid_clustered <- covid_ml %>%
  mutate(cluster = as.factor(fit$cluster))
```

**Explain:**
- `mutate()` adds a new column to the dataframe
- `as.factor()` converts cluster numbers to factor (for better plotting)
- Now we have original data + cluster assignments

### Cross-tabulation

```r
cluster_covid_table <- table(
  Cluster = covid_clustered$cluster,
  COVID_Status = covid_clustered$covid19_test_results
)
print(cluster_covid_table)
```

**Explain `table()` Function:**

**Purpose:** Creates a contingency table (cross-tabulation)

**Parameters:**
- Two or more vectors/factors
- Creates a matrix showing counts of each combination

**What the Table Shows:**
```
                COVID_Status
Cluster    Negative  Positive
      1          X         Y
      2          Z         W
```

**Interpretation:**
- If clusters align with COVID status, we'd see:
  - Cluster 1: mostly one class
  - Cluster 2: mostly the other class
- This suggests clusters found meaningful patterns!

**Teaching Points:**
- "This is validation - we're checking if unsupervised clusters match known labels"
- "Good alignment suggests the algorithm found meaningful patterns"
- "Poor alignment doesn't mean clustering failed - maybe patterns are different from COVID status"

**Metrics to Calculate (Optional):**
- **Purity:** How "pure" are clusters? (all one class = high purity)
- **Adjusted Rand Index:** Measures agreement between clusters and labels

---

## Section 7: Summary and Transition (5 minutes)

**What to Say:**
"Excellent work! Let's recap what we learned about K-means clustering."

### Key Takeaways

1. **K-means is unsupervised** - finds patterns without a target variable
2. **Choose K carefully** - use elbow method, but consider domain knowledge
3. **Scale features** - critical for distance-based algorithms
4. **Visualize clusters** - helps understand what was found
5. **Validate clusters** - check if they align with known labels (when available)

### K-means Strengths
- ✅ Simple and fast
- ✅ Works well with spherical clusters
- ✅ Easy to interpret (each point belongs to one cluster)

### K-means Limitations
- ❌ Must specify K beforehand
- ❌ Assumes clusters are spherical
- ❌ Sensitive to initialization (use high nstart!)
- ❌ Sensitive to outliers

### What's Next?

**Part 3: Random Forest**
- Supervised learning (uses target variable)
- Classification (predicts COVID test results)
- No scaling needed!
- Can handle mixed data types

**Questions to Ask:**
- "How would results change if we used K=3 or K=4?"
- "Why do you think the clusters might or might not align with COVID status?"
- "What other features might define these clusters?"

---

## Common Student Questions & Answers

**Q: How do I know if K=2 is the right choice?**
A: There's no single "right" answer. Use the elbow method as a guide, but also consider:
- Your research question (do you expect 2, 3, or more groups?)
- Interpretability (can you explain what each cluster means?)
- Cluster sizes (very small clusters might not be meaningful)

**Q: What if the elbow isn't clear?**
A: This is common! Try:
- Silhouette method (another way to choose K)
- Domain knowledge (what makes biological sense?)
- Try multiple K values and compare results

**Q: Why do we need nstart=50?**
A: K-means starts with random cluster centers. Different starts can lead to different results. Running 50 times and picking the best reduces the chance of getting a poor local optimum.

**Q: Can K-means handle non-spherical clusters?**
A: No, K-means assumes clusters are roughly spherical. For other shapes, consider:
- DBSCAN (density-based clustering)
- Hierarchical clustering
- Gaussian Mixture Models

**Q: What do the cluster centers represent?**
A: The "average" or "prototype" of each cluster. For example, if cluster 1's center has high temperature and many symptoms, that cluster represents patients with severe symptoms.

---

## Teaching Tips

1. **Emphasize the unsupervised nature:**
   - Constantly remind students we're not using the target variable
   - This is exploratory - finding hidden patterns

2. **Visualization is key:**
   - Show the elbow plot and explain how to read it
   - Use the cluster visualization to help students "see" the results

3. **Common mistakes to warn about:**
   - Forgetting to scale (show what happens if you don't!)
   - Using too low nstart (might get poor results)
   - Choosing K based only on statistics (ignore domain knowledge)

4. **Interactive elements:**
   - Have students try different K values and compare
   - Ask them to interpret what each cluster might represent
   - Compare cluster characteristics (mean values of features per cluster)

5. **Connect to biology:**
   - "What biological patterns might these clusters represent?"
   - "Do clusters align with disease severity? Risk factors?"
   - "How could this help in clinical decision-making?"
