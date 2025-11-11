#!/bin/bash

# This is a comment. Comments are ignored by the shell.

# $1 is the first argument passed to the script
echo "Counting the genes in $1"

# count the unique genes in the file
u_genes=$(gunzip -c $1 | cut -f 1 | sort -u | wc -l)

echo "There are $u_genes unique genes in $1"
