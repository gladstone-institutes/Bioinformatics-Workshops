#This script perform some basic calculations in R.
#To run this script you may select all and hit the Run button on top right of this pane ...
#... or go Cmd+A followed by Cmd+Enter on Mac. If you use Windows, ...
#... you can also go Ctrl+A followed by Ctrl+Enter.

#The following command will add 2 to 3 and store the value in variable named 'a'.
a <- sum(2,3)

#The following command will get the product of 2 and 3 and store it in 'b'.
b <- prod(2,3)

#Next, we check if a and b have equal values.
a == b

#Next, we check if a and b are not equal.
a != b
#Conclusion: Summing numbers is not the same as multiplying them!
#Time to write a paper? I think we can go to Nature or Science with this discovery.


#Check if two conditions are simultaneously true.
(a == b) & (sqrt(3) == 5)

#Check if any one of given two conditions are true.
#Vertical line is how we say 'or' in R.
(a == b) | (sqrt(3) == 5)

#Variables can also store characters.
name <- "Homo sapiens"

#Is Homo sapiens equal to human being?
name == "Human being"

#Apparently not?
