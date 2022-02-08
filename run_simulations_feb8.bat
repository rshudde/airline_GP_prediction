# first make all the directories if they don't exist
[ -d t20 ] || mkdir t20
[ -d t40 ] || mkdir t40
[ -d t60 ] || mkdir t60
[ -d t80 ] || mkdir t80

# now create the directories of the output 
[ -d RESULTS ] || mkdir RESULTS
[ -d output ] || mkdir output

# now create the data
nohup R CMD BATCH DATA_subset_datasets.R DATA_subset_datasets.out &

# now run all of the datasets for x = 20 to see how this goes yay
for i in {1..50}; do \
nohup R CMD BATCH --no-save --no-restore '--args r=i t=20' BLAH.R BLAH.out &
done