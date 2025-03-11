''' compare ice and ocn in those two folders
    ICE6G_1x1_processed: generated using C preprocess files.
    ICE6G_1x1_122_iter1: generated using Python preprocess files.
'''

import numpy as np


def compare_two_files(file1, file2):
    data1 = np.loadtxt(file1)
    data2 = np.loadtxt(file2)
    data1 = data1[:,2]
    data2 = data2[:,2]
    relative_diff = np.sqrt(np.sum((data1 - data2)**2) / np.sum(data2**2))
    return relative_diff

for epoch in range(122):
    print(f"epoch: {epoch:d}")
    # ice
    fn_ice1 = './ICE6G_1x1_processed/ice180x360.' + str(epoch)
    fn_ice2 = './ICE6G_1x1_122_iter1/ice6g180x360.' + str(epoch)
    diff_ice = compare_two_files(fn_ice1, fn_ice2)
    print(f"\t ice diff:{diff_ice:.8%}")

    # ocn
    fn_ocn1 = './ICE6G_1x1_processed/ocn180x360.' + str(epoch)
    fn_ocn2 = './ICE6G_1x1_122_iter1/ocn180x360.' + str(epoch)
    diff_ocn = compare_two_files(fn_ocn1, fn_ocn2)
    print(f"\t ocn diff:{diff_ocn:.8%}")