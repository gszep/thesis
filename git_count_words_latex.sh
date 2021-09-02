#!/bin/bash
#
# This script traverses all commits of a git repository in the current directory
# and counts the number of words that are changed, i.e. added or deleted for all
# TeX files. The output contains a time-stamp (YYYY-MM-DD). Next, there are some
# counts, viz. the number of added, deleted, and total words.
#
# You can use `gnuplot`, for example, to create nice visualizations from the raw
# data. For this, it may be useful to specify the first column as a date column:
#
#   set xdata time
#   set timefmt "%Y-%m-%d"
#
# Following this, you may plot the desired information. Use 'branch' in order to
# control for which branch the commits are being enumerated.

branch="master"
words=0

for commit in $(git rev-list --reverse $branch)
do
  date=$(git show -s --format=%cd --date=short $commit)
  added=$(git show -p --word-diff=porcelain $commit "*.tex" | grep -e '^+[^+]' | wc -w)
  deleted=$(git show -p --word-diff=porcelain $commit "*.tex" | grep -e '^-[^-]' | wc -w)
  
  words=$(($words+$added))
  words=$(($words-$deleted))
  
  echo $date $added $deleted $words
done
