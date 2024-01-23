# Exercise 1: Rscript Files -----------------------------------------------

# This is an Rscript file
# - It contains lines of code that can be executed using Ctrl+Enter
#
# - Lines that begin with # are called comments and are ignored by the compiler,
#   these can be used to document what the code does
#
# - The whole file can be executed by pressing Source above
#    *This will error unless you comment out line 49*
#
# - Refer back to the slides at the end of each Exercise

# Use Ctrl+Enter to execute each of these lines

1 + 5    # Addition

70 - 23  # Subtraction

64 / 8   # Division

8 * 8    # Multiplication

sqrt(25) # Square root

5 ^ 2    # Exponent

pi       # Built in variable for pi

# Rstudio's tab completion feature is very useful, type "sq" and press tab


# **Check your understanding**
# Calculate the area of a circle with radius = 5cm 



# End of Exercise 1 -------------------------------------------------------


# Exercise 2: Data Types and Structures -----------------------------------

##### 2.1 Data Types #####

# Different types and structures have can have different operations and
# functions done to them. Use class() to check the type of a variable

# Integer
x <- 4L
class(x)

# Numeric
y <- 1.5
class(y)

# When possible R will coerce types, for example:
class(x + y)

# It will not coerce all types though, this will throw an error
1 + "2"

# You can change the type of a variable using "as" functions
1 + as.numeric("2")

# What do you think the output of this line will be?
4 + TRUE

# What about this line?
class(c(TRUE, FALSE, TRUE, 4))

# NAs are a special case of logical, used most often to represent missing data
class(NA)

# Math operations with NAs return NA
2 + NA
mean(c(1,2,3,NA))

##### 2.2 Data Structures #####

# Defining an atomic vector
abc <- c("a", "b", "c")

# Accessing elements of a vector
abc[2]

# Range of elements
abc[1:2]

# Position of an element
which(abc == "b")

# Easily specify a vector of numbers
num <- 1:10
num

# Defining an unnamed list (note that here we overwrite the original contents of abc)
abc <- list("a", "b", 3)

# Accessing elements of unnamed list
abc[[1]]

# Defining a named list
abc <- list(a = 1, b = "2", c = 3)

# Elements of a named list can be accessed 3 ways
abc$b # Name with dollar sign
abc[["b"]] # Name with brackets
abc[[2]] # Numeric index with brackets

# Range of elements
abc[2:3]

# Numeric index to exclude elements
abc[-3]

# Defining a matrix
numeric_matrix <- matrix(1:9, nrow = 3, ncol = 3)
numeric_matrix

# Elements are accessed using row,column indexes
numeric_matrix[2, 3]
numeric_matrix[2, ] # Whole 2nd row
numeric_matrix[1:2, ] # First 2 rows

# Rows and columns of a matrix can have names
numeric_matrix <- matrix(data = 1:9,
  nrow = 3,
  dimnames = list(
    c("X", "Y", "Z"),
    c("A", "B", "C")
  )
)
numeric_matrix

# Elements can be accessed both by numeric index and column and row names
numeric_matrix[2, 3]
numeric_matrix["Y", "C"]
numeric_matrix[1:2, "C"]

# Defining a data frame
mixed_df <- data.frame(
  numbers = 1:3,
  letters = c("A", "B", "C"),
  bool = c(TRUE, FALSE, NA)
)
mixed_df

# Accessing elements of a data frame

mixed_df$numbers      # Returns the vector 
mixed_df[["numbers"]] 

mixed_df["letters"]   # Returns the column as a single column data frame

# Count rows and columns in data frames or matrices
nrow(mixed_df)
ncol(mixed_df)

# Like matrices, elements in a data frame can  also be accessed by
# numeric index or column and row names.
mixed_df[2, 3]
mixed_df["B", "letters"]
mixed_df[1:2, "letters"]

######################## **Check your understanding** ########################
## Data Types:
#
# 1) Write 2 lines of code that converts the following character vector into a
#    numeric vector and then checks the type of that vector.
vec <- c("1", "2", "3", "4", "5")


## Data Structures:
#
# 2) Write 1 line of code that accesses only elements B, C, FALSE and NA of 
#    mixed_df.
#
# Hint: line 134


# 3) Write the two ways to access only elements "2" and 3 from the list abc.
#





 



