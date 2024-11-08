if [[ -f "./main" ]]; then
  rm ./main
fi
g++ ./main.cpp -o main -O3 -fopenmp
# g++ ./main.cpp -O0 -g -o main -fopenmp
if [[ -f "./main" ]]; then
  ./main
fi
