#################################################
# Code name: a2p2_lib.py				
# Decription: Libraries for A2P2
# Developver: Kanghoon Yim (feihoom82@gmail.com)
# Date: 2020-08-01
#################################################

import sys, math, subprocess, os
from numpy import array

symbol_map = {'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, \
'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20, \
'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, \
'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, \
'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, \
'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60, \
'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, \
'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, \
'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88, 'Ac':89, 'Th':90, \
'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99, 'Fm':100, \
'Md':101, 'No':102, 'Lr':103, 'Rf':104, 'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110, \
'Rg':111, 'Cn':112, 'Uut':113, 'Uuq':114, 'Uup':115, 'Uuh':116, 'Uus':117, 'Uuo':118}
symbol_map_inv = {v: k for k, v in symbol_map.items()}
home = os.getcwd()

### Classes ###
class POSCAR(object) :
	def __init__(self,name,mag,lattice,element,natom,SD,coord,xyz,xyz_sd,xyz_tag,sym='1') :
		self.name = name
		self.mag = mag
		self.lattice = lattice
		self.element = element
		self.natom = natom
		self.SD = SD
		self.coord = coord
		self.xyz = xyz
		self.xyz_sd = xyz_sd
		self.xyz_tag = xyz_tag
		self.sym = sym

	def make_POSCAR(self) :
		poscar = ""
		poscar += self.name+'\n 1.0\n'
		for i in range(3) :
			poscar += '  %4.16F  %4.16F  %4.16F\n' %(self.lattice[i][0],self.lattice[i][1],self.lattice[i][2])
		poscar += '  '+'  '.join(self.element)+'\n'
		poscar += '  '+'  '.join([str(x) for x in self.natom])+'\n'
		if len(self.xyz_tag) != sum(self.natom) :
			self.xyz_tag = ['']* sum(self.natom)
		if self.SD == True :
			poscar += 'Selective Dynamics\n'
			poscar += self.coord+'\n'
			for i in range(len(self.xyz)) :
				poscar += '  %4.16F  %4.16F  %4.16F' %(self.xyz[i][0],self.xyz[i][1],self.xyz[i][2])
				poscar += '\t'+'   '.join(self.xyz_sd[i])+'\t'+self.xyz_tag[i]+'\n'
		else :
			poscar += self.coord+'\n'
			for i in range(len(self.xyz)) :
				poscar += '  %4.16F  %4.16F  %4.16F' %(self.xyz[i][0],self.xyz[i][1],self.xyz[i][2])
				poscar += '\t'+self.xyz_tag[i]+'\n'
		return poscar

	def chem_formula_sum(self, reduce=False, sort=False) :
		if reduce == True :
			from math import gcd
			n = self.natom[0]
			for i in range(len(self.natom)-1) :
				n = gcd(n,self.natom[i+1])
		else :
			n = 1
		cfsum = []
		for i in range(len(self.element)) :
			cfsum.append(self.element[i]+'_'+str(int(self.natom[i]/n)))
		if sort == True :
			cfsum = sorted(cfsum)
		return ' '.join(cfsum)

	def convert_coord(self, coord, recalc_xyz=False) :
		from numpy import matrix
		M1 = matrix(self.lattice)
		if self.coord[0].upper() == 'D' and coord[0].upper() == 'C' :
			# recalculate direct coordinate to remove atoms out of unit-cell 
			if recalc_xyz == True :
				for i in range(len(self.xyz)) :
					for j in range(3) :
						self.xyz[i][j] = self.xyz[i][j]%1.0
			for i in range(sum(self.natom)) :
				xyz_cvt = (matrix(self.xyz[i][0:3])*M1).getA1()
				self.xyz[i] = xyz_cvt.tolist()
			self.coord = coord
		elif self.coord[0].upper() == 'C' and coord[0].upper() == 'D' :
			for i in range(sum(self.natom)) :
				xyz_cvt = (matrix(self.xyz[i][0:3])*M1.I).getA1()
				self.xyz[i] = xyz_cvt.tolist()
			if recalc_xyz == True :
				for i in range(len(self.xyz)) :
					for j in range(3) :
						self.xyz[i][j] = self.xyz[i][j]%1.0
			self.coord = coord

	def find_sym(self, convert='original', prec=1e-3, labelling=True) :
		# convert : 'original', 'primitive', 'conventional'
		import spglib
		if self.coord[0] == 'C' :
			self.convert_coord('Direct')
		symbols = []
		for i in range(len(self.element)):
			symbol=symbol_map[self.element[i]]
			for j in range(self.natom[i]):
				symbols.append(symbol)
		cell = (self.lattice,self.xyz,symbols)
		spacegroup = spglib.get_spacegroup(cell, symprec=prec)
		print(spacegroup)
		self.sym = str(spacegroup.split('(')[1][:-1])

		if labelling == True :
			eqatoms = spglib.get_symmetry(cell)['equivalent_atoms']
			eqidx = set(eqatoms)
			label = {}
			chem = ''
			for i in eqidx :
				if chem == symbol_map_inv[symbols[i]] :
					cnt += 1
				else :
					cnt = 1
					chem = symbol_map_inv[symbols[i]]					

				label[i] = symbol_map_inv[symbols[i]]+str(cnt)
			for i in range(len(eqatoms)) :
				self.xyz_tag[i] = '!'+label[eqatoms[i]]
		
		if convert[0].upper() != 'O' and convert[0].upper() in ['P','C'] :
			if convert[0].upper() == 'P' :
				cell_sym = spglib.find_primitive(cell,symprec=prec)
			else :
				cell_sym = spglib.standardize_cell(cell,symprec=prec)
			symmetry = spglib.get_symmetry(cell_sym, symprec=prec)['equivalent_atoms']
			self.lattice = cell_sym[0].tolist()
			self.xyz = cell_sym[1].tolist()
			self.xyz_tag = []; self.xyz_sd = []
			for i in range(len(self.xyz)) :
				self.xyz_tag.append('')
				self.xyz_sd.append(['T','T','T'])
			self.coord = 'Direct'
			elem_new = []; natom_new = []
			sym_uniq = []; sym_tag = []
			for i in range(len(cell_sym[2])) :
				elem = symbol_map_inv[cell_sym[2][i]]
				if len(elem_new) == 0 or elem != elem_new[-1] :
					elem_new.append(elem)
					natom_new.append(1)
				else :
					natom_new[-1] += 1
				if symmetry[i] not in sym_uniq :
					sym_uniq.append(symmetry[i])
					cnt = 1
					for t in sym_tag :
						if elem == t[0] :
							cnt += 1
					sym_tag.append([elem,cnt])
					self.xyz_tag[i] = '!'+'_'.join([str(x) for x in sym_tag[-1]])
				else :
					self.xyz_tag[i] = '!'+'_'.join([str(x) for x in sym_tag[sym_uniq.index(symmetry[i])]])
			self.element = elem_new
			self.natom = natom_new
			if convert[0].upper() == 'C' :
				self.merge_element()

	def get_numbers(self,name) :
		numbers = []
		n = 0
		for i in range(len(self.element)) :
			if name == self.element[i] :
				for j in range(self.natom[i]) :
					numbers.append(n+j+1)
				n += self.natom[i]
			else :
				n += self.natom[i]
		return numbers

	def merge_element(self) :
		elem_m = []
		natom_m = []
		xyz_m = []
		xyz_sd_m = []
		xyz_tag_m = []
		k = 0
		for i in range(len(self.element)) :
			if not self.element[i] in elem_m :
				elem_m.append(self.element[i])
				natom_m.append(self.natom[i])
				xyz_m.append(self.xyz[k:k+self.natom[i]])
				xyz_sd_m.append(self.xyz_sd[k:k+self.natom[i]])
				xyz_tag_m.append(self.xyz_tag[k:k+self.natom[i]])
				k += self.natom[i]
			elif self.element[i] in elem_m :
				idx = elem_m.index(self.element[i])
				natom_m[idx] += self.natom[i]
				xyz_m[idx] += self.xyz[k:k+self.natom[i]]
				xyz_sd_m[idx] += self.xyz_sd[k:k+self.natom[i]]
				xyz_tag_m[idx] += self.xyz_tag[k:k+self.natom[i]]
				k += self.natom[i]
		self.element = elem_m
		self.natom = natom_m
		self.xyz = []; self.xyz_sd = []; self.xyz_tag = []
		for i in range(len(xyz_m)) :
			self.xyz += xyz_m[i]
			self.xyz_sd += xyz_sd_m[i]
			self.xyz_tag += xyz_tag_m[i]

	def atom_info(self, indices) :	# indice = [1,2,3]
		atom_info = []
		for idx in indices :
			n = 0
			for i in range(len(self.natom)) :
				n += self.natom[i]
				if idx <= n :					
					break
			elem = self.element[i]
			atom_info.append([elem,self.xyz[idx-1]])
		return(atom_info)

	def insert_atom(self, sites, verbos=True) :
	## sites = [[name,xyz],...] ex) [['Li',[0.5,0.5,0.5]]]
		for n in range(len(sites)) :
			name = sites[n][0]
			xyz = sites[n][1]
			if name in self.element :
				i = sum(self.natom[0:self.element.index(name)+1])
				self.natom[self.element.index(name)] += 1
			else :
				self.element.append(name)
				self.natom.append(1)
				i = len(self.xyz)+1
			self.xyz[i-1:i] += [xyz]
			self.xyz_sd[i-i:i] += [['T','T','T']]
			self.xyz_tag[i-1:i] += ['!'+name]
			if verbos == True :
				print(name+' atom is insert at '+str(i+1)+'th position.')

	def remove_atom(self, indices=[], condition=None, verbos=True) :
	## indices = [49,51]
	## ex) condition = 'self.xyz[2] == 0 and self.xyz[0] > 0.9'
		if condition != None :
			for i in range(len(self.xyz)) :
				if eval(condition) :
					indices.append(i+1)

		if indices == [] :
			print('You must input atomic indices or condition!')
			print("ex) remove_atom(indices=[1,2]) or remove_atom(condition='self.xyz[i][2] > 0.5')")
	
		indices = sorted(indices,reverse=True)
		i = 1
		while i < len(indices) :
			if indices[i] in indices[:i] :
				indices[i:i+1] = []
			else : i += 1;
		for idx in indices :
			n = 0
			for i in range(len(self.natom)) :
				n += self.natom[i]
				if idx <= n :					
					break
			self.xyz[idx-1:idx] = []
			self.xyz_sd[idx-1:idx] = []
			self.xyz_tag[idx-1:idx] = []
			self.natom[i] -= 1
			if verbos == True :
				print('Atom '+str(idx)+' is removed.')

		i = 0
		while i < len(self.natom) :
			if self.natom[i] == 0 :
				self.natom[i:i+1] = []
				self.element[i:i+1] = []
			else :
				i += 1

	def substitute_atom(self, chem, indices) :
		indices = sorted(indices, reverse=True)
		for idx in indices :
			subs = [[chem,self.xyz[idx-1]]]
			self.remove_atom([idx])
			self.insert_atom(subs)

	def write_atom_label(self,numbering=True) :
		n = 0
		if numbering == True :
			for i in range(len(self.element)) :
				for j in range(self.natom[i]) :
					self.xyz_tag[n] = '!'+self.element[i]+str(j+1)
					n += 1
		else :
			for i in range(len(self.element)) :
				for j in range(self.natom[i]) :
					self.xyz_tag[n] = '!'+self.element[i]
					n += 1


### Generic functions ###
def norm(x) :
	return math.sqrt(sum(i**2 for i in x))

def add_vec(vec1,vec2) :
	return([vec1[0]+vec2[0],vec1[1]+vec2[1],vec1[2]+vec2[2]])

def diff_vec(vec1,vec2) :
	vec = [vec2[x]-vec1[x] for x in range(3)]
	return([norm(vec),vec])

def diffp_vec(lattice,vec1,vec2) :
	v1 = array(vec1)
	v2 = array(vec2)
	period = array(lattice)
	min_d = 1E999
	for i in [-1,0,1] :
		for j in [-1,0,1] :
			for k in [-1,0,1] :
				vd = v2-v1+i*period[0]+j*period[1]+k*period[2]
				l = norm(vd)
				if l < min_d :
					min_d = l
					vdf = vd.tolist()[:]
	return([min_d,vdf])

def angle(a,b) :
	return math.acos(sum(a[i]*b[i] for i in [0,1,2])/(norm(a)*norm(b)))*180/math.pi

def volume(lattice) :
	a = lattice[0]; b = lattice[1]; c = lattice[2]
	return a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0])

def area(a,b) :
	return ((a[1]*b[2]-a[2]*b[1])**2+(a[0]*b[2]-a[2]*b[0])**2+(a[0]*b[1]-a[1]*b[0])**2)**0.5

def comb(nkp) :	# Make combination of path
	path=[]
	for i in range(nkp) :
		for j in range(i+1,nkp) :
			path.append([999,i+1,j+1])
	return path



### VASP functions ###
### Read structure from POSCAR file
def read_POSCAR(fposcar, input_type='file') :
	if len(fposcar.split('\n')) == 1 :
		f = open(fposcar, 'r')	
		poscar = f.read().split('\n')
		f.close()
		if len(fposcar.split('/')[-1].split('_')) > 1 :
			name = '_'.join(fposcar.split('/')[-1].split('_')[1:])
		else :
			name = 'Untitled'	

	elif len(fposcar.split('\n')) > 1 :
		poscar = fposcar.split('\n')
		name = poscar[0].strip()

	mag = float(poscar[1].strip())
	a = [float(x)*mag for x in poscar[2].split()]
	b = [float(x)*mag for x in poscar[3].split()]
	c = [float(x)*mag for x in poscar[4].split()]
	mag = 1.0
	lattice = [a,b,c]
	element = poscar[5].split()
	natom = []
	for i in range(len(element)) :
		natom[i:] = [int(poscar[6].split()[i])]
	n_tot = sum(natom)
	if poscar[7].split()[0][0] in ['S','s'] :
		SD = True
		start = 9
	else :
		SD = False
		start = 8
	if poscar[start-1].split()[0][0] in ['C','c'] :
		coord = 'Cartesian'
	elif poscar[start-1].split()[0][0] in ['D','d'] :
		coord = 'Direct'
	xyz=[]
	xyz_sd = []
	xyz_tag = []
	if SD == True :
		for i in range(start,start+n_tot) :
			xyz.append([float(x) for x in poscar[i].split()[0:3]])
			xyz_sd.append([x for x in poscar[i].split()[3:6]])
			xyz_tag.append(' '.join(poscar[i].split()[6:]))
	elif SD == False :
		for i in range(start,start+n_tot) :
			xyz.append([float(x) for x in poscar[i].split()[0:3]])
			xyz_sd.append(['T','T','T'])
			xyz_tag.append(' '.join(poscar[i].split()[3:]))
	return(POSCAR(name,mag,lattice,element,natom,SD,coord,xyz,xyz_sd,xyz_tag))

def redefine_cell(poscar, pivot='c',zaxis='c',shift=[0,0,0],recalc_xyz=False) :
	import copy
	if pivot == 'a' :
		seq = [1,2,0]
	elif pivot == 'b' :
		seq = [2,0,1]
	elif pivot == 'c' :
		seq = [0,1,2]

	coord_type = poscar.coord[0].upper()

	a = norm(poscar.lattice[seq[0]])
	b = norm(poscar.lattice[seq[1]])
	c = norm(poscar.lattice[seq[2]])
	alpha = angle(poscar.lattice[seq[1]],poscar.lattice[seq[2]])/180.*math.pi
	beta = angle(poscar.lattice[seq[2]],poscar.lattice[seq[0]])/180.*math.pi
	gamma = angle(poscar.lattice[seq[0]],poscar.lattice[seq[1]])/180.*math.pi

	cg2 = (math.cos(gamma)-math.cos(alpha)*math.cos(beta))/(math.sin(alpha)*math.sin(beta))
	sg2 = (1-cg2**2)**0.5
	cc = [0, 0, c]
	bb = [0, b*math.sin(alpha), b*math.cos(alpha)]
	aa = [a*math.sin(beta)*sg2, a*math.sin(beta)*cg2, a*math.cos(beta)]
	if pivot == 'a' :
		tmp = cc[:]
		cc = [tmp[2],tmp[0],tmp[1]]
		tmp = bb[:]
		bb = [tmp[2],tmp[0],tmp[1]]
		tmp = aa[:]
		aa = [tmp[2],tmp[0],tmp[1]]
	elif pivot == 'b' :
		tmp = cc[:]
		cc = [tmp[1],tmp[2],tmp[0]]
		tmp = bb[:]
		bb = [tmp[1],tmp[2],tmp[0]]
		tmp = aa[:]
		aa = [tmp[1],tmp[2],tmp[0]]
		
	poscar.lattice[seq[0]] = aa
	poscar.lattice[seq[1]] = bb
	poscar.lattice[seq[2]] = cc

	if zaxis == 'a' :
		nseq = [1,2,0]
	elif zaxis == 'b' :
		nseq = [2,0,1]
	elif zaxis == 'c' :
		nseq = [0,1,2]

	tmp_lattice = [[],[],[]]
	for i in range(3) :
		tmp_lattice[i] = poscar.lattice[nseq[i]]

	new_lattice = copy.deepcopy(tmp_lattice)
	for i in range(3) :
		for j in range(3) :
			new_lattice[i][j] = tmp_lattice[i][nseq[j]]

	poscar.lattice = new_lattice[:]
	for i in range(len(poscar.xyz)) :
		new_xyz = [0,0,0]
		for j in  range(3) :
			new_xyz[j] = poscar.xyz[i][nseq[j]]+shift[nseq[j]]
			if recalc_xyz == True :
				if poscar.coord[0].upper() == 'C' :
					poscar.convert_coord(coord='Direct')
					new_xyz[j] = new_xyz[j]%1.0
					poscar.convert_coord(coord='Cartesian')
				else :
					new_xyz[j] = new_xyz[j]%1.0
		poscar.xyz[i] = new_xyz[:]

### Return distance between ith and jth atom in poscar. (index start from 1)
def distance(poscar,i,j,px=True,py=True,pz=True) :
## i, j also can be position vectors
	from numpy import dot
	if poscar.coord[0] in ['D','d'] :
		lat = array(poscar.lattice)
		period = array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
	else :
		lat = array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
		period = array(poscar.lattice)

	## Check i,j whether indices or vectors
	if type(i) is int :
		pos_i = array(poscar.xyz[i-1])
	elif type(i) is list :
		pos_i = array(i)

	if type(j) is int :
		pos_j = array(poscar.xyz[j-1])
	elif type(j) is list :
		pos_j = array(j)

	vd = [1E10, 1E10, 1E10]
	vt = [0, 0, 0]
	d_min = 1E10
	if px == True :
		pxr = [-1.,0.,1.]
	else : pxr = [0]
	if py == True :
		pyr = [-1.,0.,1.]
	else : pyr = [0]
	if pz == True :
		pzr = [-1.,0.,1.]
	else : pzr = [0]
	for l in pxr :
		for m in pyr :
			for n in pzr :
				vt = pos_j-pos_i+l*period[0]+m*period[1]+n*period[2]
				d = norm(dot(vt,lat))	
				if d < d_min :
					d_min = d
					vd = vt.tolist()[:]
	return [d_min, vd]


### Analyze nearest neighbor (NN) information
def NN_info(poscar,atom_idx,LR=1.25,SR=1.05,exclude=[],tol=0.001,verbose=True,CN=True,H_on=False,vector='relative') :
# CN=True : list all the image atoms to find CN for small cell
# H_on=True : count hydrogen as nearest neighbor
	from numpy import dot
	from operator import itemgetter
	# Labeling to atoms
	if poscar.coord[0].upper() == 'D' :
		lat = array(poscar.lattice)
		period = array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
	else :
		lat = array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
		period = array(poscar.lattice)

	if type(atom_idx) is int :
		idx = atom_idx-1
		xyz = poscar.xyz[idx]
		PB = [-1.,0,1.]
	elif type(atom_idx) is list and len(atom_idx) == 3 :
		xyz = atom_idx
		PB = [-2,-1,0,1,2]
	else :
		print('Wrong parameter for atom_idx.')
		return(-1)

	LIST = []
	cnt = 0
	for i in range(len(poscar.element)) :
		for j in range(poscar.natom[i]) :
			if CN == True :
				for l in PB :
					for m in PB :
						for n in PB :
							vt = array(poscar.xyz[cnt])-array(xyz)+l*period[0]+m*period[1]+n*period[2]
							d = norm(dot(vt,lat))
							d = round(d/tol)*tol
							if vector == 'relative' :
								LIST.append([poscar.element[i],cnt+1,j+1,d,vt.tolist()])
							elif vector == 'absolute' :
								xyz_NN = [xyz[x]+vt.tolist()[x] for x in [0,1,2]]
								LIST.append([poscar.element[i],cnt+1,j+1,d,xyz_NN]) 
								# Element, atom#, elem#, distance, vector
			else :
				[d,vt] = distance(poscar,idx+1,cnt+1)
				LIST.append([poscar.element[i],cnt+1,j+1,d,vt])

			cnt += 1

	LIST = sorted(LIST, key=itemgetter(3))
	result = LIST[1:]
	k = 0
	if H_on == True :
		while result[k][0] == 'H' or result[k][3] < 0.8 or result[k][0] in exclude :
			k += 1
	else :
		while result[k][0] in exclude :
			k += 1

	min_D = result[k][3]
	cnt = 0
	NN = 0
	NN5 = 0
	NN_list = []
	while cnt < len(result) and result[cnt][3] < min_D*LR :
		if result[cnt][0] in exclude :
			cnt += 1
			continue
		if result[cnt][3] < min_D*SR :
			NN5 = NN5 + 1
		NN += 1
		NN_list.append(result[cnt])
		cnt += 1

	if verbose==True :
		if type(atom_idx) == int :
			print(f'Atom: {LIST[0][0]} ({atom_idx}) : '+', '.join([str(x) for x in poscar.xyz[atom_idx-1]]))
		else :
			print('Site: ('+', '.join([str(x) for x in atom_idx])+')')
		print(f'NN : {NN} ({LR*100}%) / {NN5} ({SR*100}%)')
		print('No.  Element\tAtom #\tElem #\tDistance#\tVector')
		for i in range(len(NN_list)) :
			print(f'{i+1:3d}  '+'\t'.join([str(x) for x in NN_list[i]]))

	return NN_list
