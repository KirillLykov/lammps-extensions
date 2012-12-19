#! /bin/bash
# requirements - Velocities or Bonds section is after Atoms
# Dihedrals Section is right after Angles

start_atoms=`grep -n "Atoms" $1 | cut -d ":" -f 1`
start_atoms=$(( $start_atoms + 2))
end_atoms=`grep -n "Velocities" $1 | cut -d ":" -f 1`

if ! [[ "$end_atoms" =~ ^[0-9]+$ ]] ; then
   end_atoms=`grep -n "Bonds" $1 | cut -d ":" -f 1`
fi
end_atoms=$(($end_atoms - 2))

start_angles=`grep -n "Angles" $1 | cut -d ":" -f 1`
start_angles=$(( $start_angles + 2))
end_angles=`grep -n "Dihedrals" $1 | cut -d ":" -f 1`
end_angles=$(($end_angles - 2))
#$(cat $1 | wc -l) - get length of file

atoms_count=$(($end_atoms - $start_atoms + 1))
angles_count=$(($end_angles - $start_angles + 1))

filename=$(basename "$1")
extension="${filename##*.}"
filename="${filename%.*}"
outfilename=$filename".plt"

header="TITLE     = \"$filename\"
VARIABLES = \"X\"
\"Y\"
\"Z\"
ZONE T=\"ZONE 001\"
 STRANDID=0, SOLUTIONTIME=0
  Nodes=$atoms_count, Elements=$angles_count, ZONETYPE=FETriangle
   DATAPACKING=POINT
    DT=(SINGLE SINGLE SINGLE )"

echo "$header" > $outfilename
sed -n "$start_atoms,$end_atoms""p" $1 | sort -k1 -n | cut -d " " -f 4-6 >> $outfilename
sed -n "$start_angles,$end_angles""p" $1 | cut -d " " -f 3-5 >> $outfilename

