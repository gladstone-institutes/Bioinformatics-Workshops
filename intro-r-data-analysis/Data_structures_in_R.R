#Vectors.
#This is a numeric vector.
a <- c(3, 5, 2.5, 0, 9)
mean(a)

#Vectors (as well as other data objects) may contain NA.
#NA stands for Not Available.
a <- c(3, 5, 6.8, NA)
mean(a)

#To ignore NA while taking mean.
mean(a, na.rm = TRUE)

#Vectors are ordered sequences.
#This is a character vector.
b <- c("Homo sapiens", "Martians", "Blue.Whales", "Homo sapiens")
table(b)

#----------
#Data frame
#Tabular data possibly with mix of numeric, character, factor or logical type entries.
#Example: Iris setosa data that we worked with.
df <- data.frame(a = 1:4,
                 name = b)

#---------
#Matrix
#Matrices are tabular like data frames but store only one type of entries.
df_matrix <- as.matrix(df)

#May convert between data structures.
df_numeric <- as.numeric(df_matrix)
df_character <- as.character(df_numeric)

#----------
#List
#Lists are flexible data structures.
#Lists have named fields which can contain arbitrary data—vectors, other lists, strings, functions, and anything else.
#We may wish to keep related data together in one place. 
#Say you search for a sequence in genome using an R library.
#Your ideal output might have mixed entries. For example,
#if sequence is found on a chromosome, you want table containing start and end loci on the 
#chromosome and matched length.
#If not found, you want a string that says "Not found."
result <- list(chr1 = data.frame(Start = c(5, 100, 200),
                                 End = c(70, 150, 230),
                                 Length = c(66, 50, 30)
                                 ),
               chr2 = "Not found."
               )