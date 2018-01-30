# This script moves all PDF files to a folder called pdf
# Create the pdf folder if you haven't already
BASEDIR=$(dirname $0)
echo "Script location: ${BASEDIR}"
cd ${BASEDIR}
mv *.pdf pdf
exit 0