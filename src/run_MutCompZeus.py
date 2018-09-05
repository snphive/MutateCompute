from multiprocessing.dummy import Pool as ThreadPool
import sys
import time
import multiprocessing as mp
from src.MutateFasta import MutateFasta


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
