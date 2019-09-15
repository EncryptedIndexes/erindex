# Change the following parameters at your convenience
# Absolute path of the directory containing the ERIndex executable
ERINDEX_BIN_DIR=/home/fernando

# Absolute path of your database root (create it before launching this script)
DATABASE_ROOT=/home/fernando/testdb

$ERINDEX_BIN_DIR/ERIndex buildindex $DATABASE_ROOT ../data/sequences/chr20_1MB_rand_coll.fa 20_1MB_rand.eri 20 
