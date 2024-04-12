# python3 experiment/benchmark.py 3 ~/graph-data/ncuts-graphs/

# python3 experiment/comparison.py ~/graph-data/chaco/ncuts-graphs/

# python3 experiment/benchmark.py 3 ~/graph-data/our-test-suite/
python3 experiment/threshold.py 3 ~/graph-data/chaco/ncuts-graphs/ -c -n 10
python3 experiment/threshold.py 3 ~/graph-data/chaco/duplicates/ -c -n 10
