from multiprocessing.dummy import Pool as ThreadPool
import sys
import time
import multiprocessing as mp
from src.MutateFasta import MutateFasta

__author__ = "Shahin Zibaee"
__copyright__ = "Copyright 2018, The Switch lab, KU Leuven"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Shahin Zibaee"
__email__ = "shahinzibaee@hotmail.com"
__status__ = "Development"


def main(arg_list=sys.argv[1]):
    print('This computer has ' + str(mp.cpu_count()) + ' cpus')
    start_time = time.perf_counter()
    pool = ThreadPool()
    results = pool.map(func=MutateFasta.mutate_every_residue_to_every_aa_write_1_file, iterable=arg_list)
    pool.close()
    pool.join()
    elapsed = time.perf_counter() - start_time
    print('Time taken to complete iterating through fasta_list (after methods have been called): ' + str(elapsed))
    return results


if __name__ == '__main__':
    main()
