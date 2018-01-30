# This script moves all PDF files to a folder called pdf
# Create the pdf folder if you haven't already
MAINLATEX="main"
TODAY=`date +%Y-%m-%d`
BASEDIR=$(dirname $0)
cd ${BASEDIR}
latexdiff --exclude-textcmd="section,subsection" $MAINLATEX\_flat\_old.tex $MAINLATEX\_flat.tex > $MAINLATEX\_diff\_$TODAY.tex
exit 0