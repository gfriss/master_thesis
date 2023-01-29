cd $1

for sim in *; do
	if [ -d "$sim" ]; then
		cd $sim
		rm *out
		cd ..
	fi
done
