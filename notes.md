
# exact file names that depend on other software, might be error prone in case
# files are named wrong

grep -nr "import " | grep -v '%' | grep -v '=' | cut -f 3 -d ':' | cut -f 1 -d ';' | sed 's/\<import\>//g' | tr -d '  ' |sed 's/^/+/' | sed 's:\.:\/+:g' |

# This copies all the required files from the big repository
'''bash
cat files.txt | while read line; do \
mkdir -p bin/$line \
directory='somepathonyourpc'
cp -R 'directory$line'.m' \
bin/$line'.m'; done
'''



cd src/;grep -nr "import " | grep -v '%' | grep -v '=' | cut -f 3 -d ':' | cut -f 1 -d ';' | sed 's/\<import\>//g' | tr -d [:blank:]' ' |sed 's/^/+/' | sed 's:\.:\/+:g' | grep -vw "+Hca"  | grep -v "'" | sort | uniq | sed 's/\(.*\)+/\1/' > files.txt


cd ../; cat src/files.txt | while read line; do line2="${line%/*}/";mkdir -p bin/$line2; cp -nR '/home/albyback/git/Projects/hca/hca_tests/hca_3.2/lldev/src/MATLAB/'$line'.m' bin/$line'.m' || true; done
rm src/files.txt;

cd bin/;
# should loop this until binfiles does not change
for i in {1..10}
do
	grep -nr "import " | grep -v '%' | grep -v '=' | cut -f 3 -d ':' | cut -f 1 -d ';' | sed 's/\<import\>//g' | tr -d [:blank:]' '  |sed 's/^/+/' | sed 's:\.:\/+:g' | grep -vw "+Hca"  | grep -v "'" | sort | uniq | sed 's/\(.*\)+/\1/'  > binfiles.txt
 	cat binfiles.txt | while read line; do line2="${line%/*}/";mkdir -p $line2; cp -nR '/home/albyback/git/Projects/hca/hca_tests/hca_3.2/lldev/src/MATLAB/'$line'.m' $line'.m' || true; done
done


Export information scores into a txt file. Ask Ville for examples
