# Change the following parameters at your convenience
# Absolute path of the directory containing the ERIndex executable
ERINDEX_BIN_DIR=/home/fernando

# Absolute path of your database root (create it before launching this script)
DATABASE_ROOT=/home/fernando/testdb

$ERINDEX_BIN_DIR/ERIndex addref $DATABASE_ROOT 20 ../data/references/hs37d5_chr20.fa 
