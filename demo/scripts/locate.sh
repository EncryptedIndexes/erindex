# Change the following parameters at your convenience
# Absolute path of the directory containing the ERIndex executable
ERINDEX_BIN_DIR=/home/fernando

#Before launching this command you must change some parameters in the XML file (i.e. the database root)
$ERINDEX_BIN_DIR/ERIndex locate ../patterns/20_1MB_rand.xml > ../results/locate_20_1MB_rand_result.txt
