# This script calls the Python script flatex.
# Make sure flatex is already installed.
# Edit this script to rename what the main LaTeX file is called.
MAINLATEX="main"  # without extension
BASEDIR=$(dirname $0)
cd ${BASEDIR}
flatex $MAINLATEX.tex $MAINLATEX\_flattened.tex
exit 0