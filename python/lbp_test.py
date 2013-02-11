import cv2
import numpy as np
import sys

sys.path.insert(0, '/home/pmor6790/Code/External/LBP/build/python')

import lbp


if __name__ == '__main__':

    img = cv2.imread('../../test_image_1.bmp', 0)
    img = img.astype(np.float64)

    l = lbp.LBP(8, lbp.LBP_MAPPING_HF)
    l.calcLBP(img)

    hist = l.calcHist().getHist()

    print(hist)
