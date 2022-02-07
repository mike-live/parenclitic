from IG_split import IG_split, IG_split_old
import sys
import numpy as np
from scipy import stats 
from kde import gaussian_kde as gaussian_kde_new
from my_cov import my_cov

n = 87
np.random.seed(10)
p = np.random.rand(n)
y = np.random.randint(2, size=(n,)) * 2 - 1

X_i = [0.0468271 , 0.0504004 , 0.06022958, 0.04297324, 0.04669779,
       0.0498247 , 0.06717592, 0.04643571, 0.0516535 , 0.04146632,
       0.05070179, 0.0541143 , 0.0456006 , 0.05096081, 0.05963892,
       0.03480463, 0.04535078, 0.05626587, 0.03779833, 0.04683609,
       0.05646019, 0.03229787, 0.04422927, 0.05926363, 0.04151269,
       0.04886416, 0.04882487, 0.04909128, 0.05308635, 0.06781207,
       0.04943665, 0.0652729 , 0.04567274, 0.04028068, 0.06629099,
       0.05823978, 0.0583286 , 0.05544021, 0.04717778, 0.04713915,
       0.04889964, 0.04736223, 0.05462356, 0.06185693, 0.03658799,
       0.04443223, 0.05801454, 0.03885529, 0.04412721, 0.03724905,
       0.03190382, 0.03093849, 0.04662203, 0.04212832, 0.05009559,
       0.03863529, 0.03645103, 0.04402535, 0.05372237, 0.04914704,
       0.05634273, 0.05524715, 0.04201493, 0.05984507, 0.04574849,
       0.06333333, 0.04344941, 0.05422289, 0.04250352, 0.03983114,
       0.04195737, 0.0582264 , 0.05297043, 0.05578824, 0.04538403,
       0.0461121 , 0.05496988, 0.04610783, 0.04088443, 0.05425539,
       0.03905082, 0.04544259, 0.04602544, 0.0382041 , 0.04828045,
       0.05205117, 0.04216909]
	   
X_j = [0.01547616, 0.01590072, 0.02072823, 0.02125974, 0.01667862,
       0.01770266, 0.02047575, 0.01679467, 0.01608288, 0.01566297,
       0.0182412 , 0.01605436, 0.01493159, 0.01731028, 0.0160624 ,
       0.01612139, 0.01626472, 0.01592707, 0.01650127, 0.01752987,
       0.01603506, 0.01591711, 0.01460928, 0.0160055 , 0.0157178 ,
       0.01481787, 0.01471465, 0.01669729, 0.01685874, 0.01805606,
       0.01590859, 0.02168766, 0.02422972, 0.01690858, 0.01714365,
       0.01771006, 0.01581072, 0.01662557, 0.01524931, 0.01504924,
       0.01567862, 0.01455504, 0.01563879, 0.01600705, 0.01496816,
       0.01603496, 0.01608734, 0.01545972, 0.01800516, 0.01641634,
       0.01622168, 0.0173019 , 0.01594049, 0.0153755 , 0.01767455,
       0.01555547, 0.01578296, 0.01510381, 0.01663375, 0.01759583,
       0.01693599, 0.02079662, 0.02102147, 0.0189431 , 0.01759507,
       0.01482431, 0.0171931 , 0.015842  , 0.0266297 , 0.01549387,
       0.01566862, 0.01510952, 0.01444809, 0.01721404, 0.01650868,
       0.01826619, 0.01570751, 0.01678429, 0.01441635, 0.0160326 ,
       0.01597796, 0.01659173, 0.01803947, 0.01740419, 0.0166103 ,
       0.01588576, 0.01746806]	   

X_i = np.array(X_i)
X_j = np.array(X_j)

y[:29] = 1
y[29:29 * 2] = -1
y[29 * 2:29 * 3] = -2

mask = y.copy()
#print(mask.shape, mask.dtype)
#print(X_i.shape, X_i.dtype)
#print(X_j.shape, X_j.dtype)


X_prob_i, X_prob_j = X_i[mask == -1], X_j[mask == -1]
data = np.array([X_prob_i, X_prob_j])
kde = gaussian_kde_new(data)

print(my_cov(data))

data = np.array([X_i, X_j])
p = np.array(kde(data))
fit_mask = (mask == -1) | (mask == 1)

r1 = IG_split(p[fit_mask], mask[fit_mask])
r2 = IG_split_old(p[fit_mask], mask[fit_mask])
print(r1, r2)

X_prob_i, X_prob_j = X_i[mask == -1], X_j[mask == -1]
data = np.array([X_prob_i, X_prob_j])
kde = stats.gaussian_kde(data)

print(np.cov(data))

data = np.array([X_i, X_j])
p = np.array(kde(data))
fit_mask = (mask == -1) | (mask == 1)

r1 = IG_split(p[fit_mask], mask[fit_mask])
r2 = IG_split_old(p[fit_mask], mask[fit_mask])
print(r1, r2)


sys.exit(0)


import numpy as np

data = np.random.rand(2, 100)

import timeit
a = np.cov(data)
print(a)

from numba import float64, int32, int64, uint64, int8, jit

class numba_config:
    cache = False
    nopython = True
    nogil = True
    parallel = True

from my_cov import my_corrcoef, my_cov

a = my_cov(data)
#print(a)

print(my_corrcoef(data))
print(np.corrcoef(data))


np.cov = my_cov

#print(timeit.timeit("a = np.cov(data)", setup="from __main__ import np, data", number=20000))
#print(timeit.timeit("a = my_cov(data)", setup="from __main__ import np, data, my_cov", number=400000))


a = np.cov(data)
print(a)


'''from kde import gaussian_kde
import numpy as np
data = np.array([[1, 2, 3, 4], [3, 5, 21, 8]], dtype=np.float32)
kde = gaussian_kde(data)

data = np.array([[1, 2, 3, 4, 5, 6], [3, 5, 21, 8, 7, 8]], dtype=np.float32)
res = kde(data)

print(np.mean(res))

'''