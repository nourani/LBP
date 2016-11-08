import cv2
import numpy as np
import os
import lbp


if __name__ == '__main__':

    filename = os.path.join(os.path.dirname(__file__), '../test_image_1.pgm')

    img = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
    #img = img.astype(np.float64)

    l = lbp.LBP(8, lbp.LBP_MAPPING_HF)
    l.calcLBP(img, 1, True) # We need to set borderCopy to true to keep the size of the image

    labelMask = np.zeros(img.shape[:2], dtype="uint8")
    labelMask[1:20,1:20] = 255
    print 'image shape=', img.shape
    print 'mask shape =', labelMask.shape
    
    hist = l.calcHist(labelMask).getHist()

    print 'hist=', hist
    #print(l.getLBPImage())
