### Octa_info v1.0 ###
### Developers: feihoom82@gmail.com, kyjung1020@gmail.com, tndla0629@gmail.com
### Last update : 2024-09-03

from lib.base import POSCAR, read_POSCAR, NN_info, diff_vec, distance
from lib.cif2vasp import cif2vasp
import sys, math, argparse
import numpy as np

# Check octahedral and their connections
def chk_octa(P, B, X, LR=1.35, verbose=True) :
	B_idx = []
	X_idx = []
	for bb in B :
		B_idx += P.get_numbers(bb)
	for xx in X :
		X_idx += P.get_numbers(xx)

	excl = []; excl2 = []
	for elem in P.element :
		if elem not in B :
			excl.append(elem)	# exclude except B

	for elem in P.element :
		if elem not in X :
			excl2.append(elem)	# exclude except X

	B_done = []
	X_none = []; X_corner = []; X_edge = []; X_face = []; X_unknown = []

	X_Ball = []

	tot = 0
	for i in B_idx :
		# Check every nearest B - B combination
		B2B = NN_info(P, i, exclude=excl, LR=1.75, vector='absolute', verbose=False)
		# Check NN of Bi
		Bi = NN_info(P, i, exclude=excl2, LR=LR, vector='absolute', verbose=False)

		CN_Bi = len(Bi)

		if CN_Bi != 6 :
			print(f'Not all {B}-{X} bonds are octahedral!')
			return({'none':'null','corner':'null','edge':'null','face':'null','unclassified':'null'})

		# Check B-X are in octahedral connection
		for X1 in Bi :
			if X1[1] not in X_Ball :	# Count all X coordinating B
				X_Ball.append(X1[1])
			X2X = []
			for X2 in Bi :
				if X1[4] != X2[4] != 0 :
					X2X.append(distance(P,X1[4],X2[4],px=False,py=False,pz=False)[0])
			X2X = sorted(X2X)
			cnt = 0
			for j in X2X :
				if j < X2X[0]*(2**0.5):
					cnt += 1
			if cnt < 4 :
				print(f'{B}({i}th)-{X} bonds are not octahedral! (Still, CN=6)')
				return({'none':'null','corner':'null','edge':'null','face':'null','unclassified':'null'})

		### Filter distorted structure ###
		for j in range(len(B2B)) :	# For every nearest B from Bi (Bj)
			if B2B[j][1] in B_done :
				continue
			Bj_img = B2B[j][4]
			# Check NN of Bj
			Bj = NN_info(P, Bj_img, exclude=excl2, LR=LR, vector='absolute', verbose=False)
			cnt = 0
			X_share = []    # Reset counting shared X 
			# Count same X site that Bi and Bj share
			for n in range(len(Bi)) :
				for m in range(len(Bj)) :
					i_img = [round(x,6) for x in Bi[n][4]]
					j_img = [round(x,6) for x in Bj[m][4]]
					if i_img == j_img :
						cnt += 1
						if Bi[n][1] not in X_share :	# Count shared X
							X_share.append(Bi[n][1])

			if cnt == 1 :
				for t in X_share :
					if not t in X_corner :
						X_corner.append(t)
			elif cnt == 2 :
				for t in X_share :
					if not t in X_edge :
						X_edge.append(t)
			elif cnt == 3 :
				for t in X_share :
					if not t in X_face :
						X_face.append(t)
			elif cnt > 3 :
				for t in X_share :
					if not t in X_unknown :
						X_unknown.append(t)

		B_done.append(i)

	for i in X_Ball :
		if not i in X_corner+X_edge+X_face+X_unknown :
			X_none.append(i)

	w_none = len(X_none)/len(X_Ball)
	w_corner = len(X_corner)/len(X_Ball)
	w_edge = len(X_edge)/len(X_Ball)
	w_face = len(X_face)/len(X_Ball)
	w_unknown = len(X_unknown)/len(X_Ball)


	X_other = []
	for i in X_idx :
		if not i in X_Ball :
			X_other.append(i)

	if verbose == True :
		print('Octa-connection :')
		print('\tNone\tCorner\tEdge\tFace\tunclassified')
		print(f'count\t{len(X_none)}\t{len(X_corner)}\t{len(X_edge)}\t{len(X_face)}\t{len(X_unknown)}')
		print(f'(%)\t{round(w_none*100)}%\t{round(w_corner*100)}%\t{round(w_edge*100)}%\t{round(w_face*100)}%\t{round(w_unknown*100)}%')
	if len(X_other) != 0 :
		print(f'(*{round(len(X_other)/len(X_idx)*100)}% X anions are not coordinating B.)')

	#return({'none':X_none,'corner':X_corner,'edge':X_edge,'face':X_face,'unclassified':X_unknown})
	return({'none':w_none,'corner':w_corner,'edge':w_edge,'face':w_face,'unclassified':w_unknown})


# 가장 멀리 떨어진 점 찾기
def farthest_point(points, base_point) :
	max_distance = 0
	farthest_point = None
	for point in points :
		if point == base_point :
			continue
		d = diff_vec(base_point, point)[0]
		if d > max_distance:
			max_distance = d
			farthest_point = point
	return farthest_point


# 벡터만들고 사면체 부피 계산
def volume(A, B, C, D) :
	AB = np.array([B[0] - A[0], B[1] - A[1], B[2] - A[2]])
	AC = np.array([C[0] - A[0], C[1] - A[1], C[2] - A[2]])
	AD = np.array([D[0] - A[0], D[1] - A[1], D[2] - A[2]])

	cross_product = np.cross(AB, AC)
	dot_product = np.dot(cross_product, AD)

	volume = abs(dot_product*(1/6))
	return volume


# B - X Octa의 Distortion vector index 계산
def distortion_vec_index(P, B, X) :
	P.convert_coord(coord='Cartesian')

	B_idx = []
	for bb in B :
		B_idx += P.get_numbers(bb)

	values = [] # index values for each octahedral

	excl2=[] #revised
	for elem in P.element : #revised
		if elem not in X : #revised
			excl2.append(elem) #revised

	for i in B_idx :
		nn = NN_info(P, i, LR=1.35, exclude=excl2, verbose=False) #revised

		if len(nn) == 6 :
			Sum = [0, 0, 0]
			length = 0
			for j in range(len(nn)) :
				for k in [0,1,2] :
					Sum[k] += nn[j][4][k]
				vec_norm = np.linalg.norm(Sum)
				length += nn[j][3]
			values.append(vec_norm/(length/(len(nn))))

	if len(values) == len(B_idx) :
		Dist_vec_index = sum(values)/len(values)
	else :
		Dist_vec_index = 'Incomplete'

	return [Dist_vec_index, values]
	

# B - X Octa의 Quadratic elongation 계산
def quadratic_elongation(P, B, X) :
	P.convert_coord(coord='Cartesian')
	B_idx = []
	for bb in B :
		B_idx += P.get_numbers(bb)

	values = []


	excl2=[] #revised
	for elem in P.element : #revised
		if elem not in X : #revised
			excl2.append(elem) #revised

	for i in B_idx :
		LIST = []
		nn = NN_info(P, i, LR=1.35, exclude=excl2, verbose=False) #revised
		if len(nn) == 6:
			for j in range(len(nn)) :
				LIST.append(nn[j][4])

			#첫 번째 좌표에서 가장 먼 점을 찾고 리스트에서 제거
			FP=farthest_point(LIST,LIST[0])
			LIST.remove(FP)

			#나머지 좌표 중 두 번째 좌표와 가장 먼 점 찾기
			FP2=farthest_point(LIST,LIST[1])
			LIST.remove(FP2)

			V1 = volume(LIST[0], LIST[1], LIST[2], FP2)
			V2 = volume(LIST[0], LIST[1], LIST[3], FP2)
			V3 = volume(FP, LIST[1], LIST[2], FP2)
			V4 = volume(FP, LIST[1], LIST[3], FP2)
			total_V = V1+V2+V3+V4

			### Quadratic Elongation 구하기 ###
			l0 = (total_V*(3/4)) ** (1/3)
			total_K = 0
			for k in range(len(nn)) :
				K = (nn[k][3]/l0) ** 2
				total_K += K
				QE = (1/6)*total_K
			values.append(QE)


	if len(values) == len(B_idx) :
		sum_QE = sum(values)
		avg_QE = sum_QE/len(values)  #revised
	else :
		avg_QE = 'Incomplete'

	return [avg_QE, values]

### Main ###
if __name__ == '__main__' :
	parser = argparse.ArgumentParser()
	parser.add_argument("filename")
	parser.add_argument("elem_comb", nargs='*', help="[element1,element2,...] - [element3,element4,...]")
	parser.add_argument("-c", "--cutoff_ratio", type=float, default=1.35, help="Cut-off ratio for counting nearest neighbors based on the shortest bond. (default=1.35)")
	args = parser.parse_args()

	try :
		if args.filename.split('.')[-1].upper() == 'CIF' :
			P = cif2vasp(args.filename, primitive=False)
			P.merge_element()
		else :
			P = read_POSCAR(args.filename)

		BX = (''.join(args.elem_comb)).replace('[','').replace(']','')
		B = BX.split('-')[0].split(',')
		X = BX.split('-')[1].split(',')
	except :
		print("Wrong input expression.\nSee help : python", sys.argv[0], "--help")
		sys.exit(-1)

	print(f'{B}-{X}')
	
	LR = args.cutoff_ratio
	ret = chk_octa(P, B, X, LR=LR)

	if ret != -1 :
		Dist_vec_index = distortion_vec_index(P, B, X)[0]
		avg_QE = quadratic_elongation(P, B, X)[0]
	else :
		Dist_vec_index = '-'
		avg_QE = '-'

	print(f'\nDistortion vector index (avg) : {Dist_vec_index}\nQuadratic Elongation (avg) : {avg_QE}\n')
