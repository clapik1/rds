BINARY="/home/clapik/workspace/rds/bin/Debug/rds"
MESHES="/home/clapik/workspace/meshes/*.msh2"
OUTPUT_FOLDER="/home/clapik/workspace/temp"

for MESH in $MESHES
do
    NAME=$(basename $MESH | cut -d'.' -f 1)
    $BINARY $MESH $OUTPUT_FOLDER/$NAME.dat
done
