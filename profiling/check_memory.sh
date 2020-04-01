#/bin/sh

ifile="../test/regression/models.txt"

while IFS= read -r line
do
      echo ""
      echo ""
      echo "================================================================================"
      echo ""
      echo "Beginning test: $line"
      echo ""
      echo ""
      valgrind --leak-check=full --track-origins=yes ../util/cxx_interface/cxxsimple ../test/regression/reference.xml "$line" 0.02 100.0 20 300
      echo ""
      echo ""
      echo "================================================================================"
done < "$ifile"

