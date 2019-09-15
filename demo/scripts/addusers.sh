# Change the following parameters at your convenience
# Absolute path of the directory containing the ERIndex executable
ERINDEX_BIN_DIR=/home/fernando

# Absolute path of your database root (create it before launching this script)
DATABASE_ROOT=/home/fernando/testdb

$ERINDEX_BIN_DIR/ERIndex adduser $DATABASE_ROOT user_1 
$ERINDEX_BIN_DIR/ERIndex adduser $DATABASE_ROOT user_2 
$ERINDEX_BIN_DIR/ERIndex adduser $DATABASE_ROOT user_3 
