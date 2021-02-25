import cv2
from skimage.metrics import structural_similarity
from skimage.metrics import mean_squared_error

# Read images from file.
im1 = cv2.imread('./true.png')
im2 = cv2.imread('./test.png')

#tempDiff = cv2.subtract(im1, im2)

gray1 = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
gray2 = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)

mse = mean_squared_error(gray1, gray2)
score = structural_similarity(gray1, gray2)

print(f"MSE : {mse: .5f}")
print(f"SSIM: {score: .5f}")