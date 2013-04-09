import cv2
import numpy as np
import os
import lbp


if __name__ == '__main__':

    filename = os.path.join(os.path.dirname(__file__), '../test_image_1.bmp')

    img = cv2.imread(filename, 0)
    img = img.astype(np.float64)

    l = lbp.LBP(8, lbp.LBP_MAPPING_HF)
    l.calcLBP(img)

    hist = l.calcHist().getHist()

    print(hist)
    print(l.getLBPImage())
