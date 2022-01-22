import numpy as np
import pandas as pd
import random
import constants
import os
from multiprocessing import Pool, Manager
import argparse


def permute_slice(original, slice_size=2, seed=1):
    
    local_r_state = np.random.RandomState(seed)
    r_rows = local_r_state.randint(0, original.shape[0], slice_size)
    r_cols = local_r_state.randint(0, original.shape[1], 1) 
    original[r_rows, r_cols] = original[np.flip(r_rows), r_cols]
         
    return original


def calc_a_scores(mat):
    mat = np.array(mat)
    mat[mat == -1] = 0
    return np.sum(mat, axis=1)


def calc_d_scores(mat):
    mat = np.array(mat)
    mat[mat == 1] = 0
    return np.abs(np.sum(mat, axis=1))


def calc_e_score(mat, o_a_scores, o_d_scores, n_iteration):
    cur_a_score = calc_a_scores(mat)
    cur_d_score = calc_d_scores(mat)

    t_a = 0.002 + 0.0001 * (n_iteration / 10000)  # (n_iteration/100+1)/100.0
    t_d = 0.002 + 0.0001 * (n_iteration / 10000)  # (n_iteration/100+1)/100.0

    exp_a = np.sum((cur_a_score - o_a_scores) / (o_a_scores + 1.0))
    exp_d = np.sum((cur_d_score - o_d_scores) / (o_d_scores + 1.0))
    # if n_iteration % 100000 == 0:
    #     print "iteration: ", n_iteration, " sa expressions: ", exp_a + exp_d
    return t_a * exp_a + t_d * exp_d


def permute_matrix(args):

    original, file_name, sa_factor, n_iteration, arr = args    
    original = np.array(original)
    original_c = np.array(original)

    o_a_scores = calc_a_scores(original_c)
    o_d_scores = calc_d_scores(original_c)

    for a in range(sa_factor / 10):
        original = permute_slice(original, seed=a+(sa_factor / 10)*n_iteration)

    prev_e_score = -1
    cur_e_score = None
    
    for a in range(sa_factor):
        original_c = permute_slice(np.array(original), seed=a+sa_factor)       
        e_score = calc_e_score(original_c, o_a_scores, o_d_scores, a + 1)
        # if a % 100000 == 0:  # or ( a> 5000 and 1-e_score > 0):
        #     print "{}	chance for accepting the step: {}".format(n_iteration, 1 - e_score)
        np.random.seed(int(random.random()*10000))
        step = np.random.binomial(1, np.min([1.0, np.max([0.0, 1 - e_score])]))

        if step: #  or cur_e_score < e_score: # e_score < prev_e_score:
            original = original_c
            cur_e_score = e_score
        if a % 100000 == 0 and a!=0:
          print "file: {} n_matrix: {} cur_e_score: {}".format(file_name, n_iteration, cur_e_score)
        prev_e_score = e_score

    # print "permuted", np.sum(original == original_c), "total: {}".format(original.size)

    try:
        os.makedirs(os.path.join(constants.CACHE_FOLDER, file_name))
    except:
        pass

    np.save(
        open(os.path.join(constants.CACHE_FOLDER, file_name, "{}_perm_{}.tsv".format(file_name, n_iteration)), 'w+'),
        original)
    try:
       arr.append(1)
       print "done {} permutations".format(np.sum(arr))
    except:
       pass
    return original


def generate_permuted_matrices(file_name, n_start, n_end, p_factor, sa_factor, use_cache):
    mat = pd.read_csv(os.path.join(constants.DATASETS_FOLDER, "{}.tsv".format(file_name)), sep='\t', index_col=0)
    p = Pool(p_factor)
    arr=Manager().list([])
    params = []
    mat[pd.isna(mat)] = 0  # MODIFY THIS LINE
    for a in np.arange(n_start,n_end):
        if use_cache and os.path.exists(os.path.join(constants.CACHE_FOLDER, file_name, "{}_perm_{}.tsv".format(file_name, a))):
            arr.append(1)
        else:
            params.append([mat, file_name, sa_factor, a, arr])

    print "permuting {}/{} matrices ({} exist in cache)".format(n_end-n_start-len(arr), n_end-n_start, len(arr))
    p.map(permute_matrix, params)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--p_factor', dest='p_factor', default="8")
    parser.add_argument('--sa_factor', dest='sa_factor', default="10000")
    parser.add_argument('--n_start', dest='n_start', default="0")
    parser.add_argument('--n_end', dest='n_end', default="20")
    parser.add_argument('--file_names', dest='file_names', default="brca_wgd_plus")
    parser.add_argument('--use_cache', dest='use_cache', default="true")
  
    args = parser.parse_args()

    p_factor = int(args.p_factor)
    sa_factor = int(args.sa_factor)
    n_start= int(args.n_start)
    n_end= int(args.n_end)
    file_names=args.file_names.split(",")
    use_cache=args.use_cache=="true"

    for cur_file in file_names:
        generate_permuted_matrices(cur_file, n_start, n_end, p_factor, sa_factor, use_cache)
    


