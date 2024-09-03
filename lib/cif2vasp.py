# AUTHOR: Kanghoon Yim (khyim@kier.re.kr)
# FILE:cif2vasp.py
# Version:1.0
# DATE:2019-5-14

import sys, math, cvt_cell
from operator import itemgetter, attrgetter
from lib.base import POSCAR, read_POSCAR

def cif2vasp(cif, primitive=False) :
	### Atomic symbols <-> numbers
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

	### Symmetry k-points table number ###
	k_table = {'1':'19', '2':'19', '3':'15', '4':'15', '5':'16', '6':'15', '7':'15', '8':'16', '9':'16', '10':'15', \
'11':'15', '12':'16', '13':'15', '14':'15', '15':'16', '16':'7', '17':'7', '18':'7', '19':'7', '20':'11', \
'21':'11', '22':'8', '23':'10', '24':'10', '25':'7', '26':'7', '27':'7', '28':'7', '29':'7', '30':'7', \
'31':'7', '32':'7', '33':'7', '34':'7', '35':'11', '36':'11', '37':'11', '38':'11', '39':'11', '40':'11', \
'41':'11', '42':'9', '43':'8', '44':'10', '45':'10', '46':'10', '47':'7', '48':'7', '49':'7', '50':'7', \
'51':'7', '52':'7', '53':'7', '54':'7', '55':'7', '56':'7', '57':'7', '58':'7', '59':'7', '60':'7', \
'61':'7', '62':'7', '63':'11', '64':'11', '65':'11', '66':'11', '67':'11', '68':'11', '69':'9', '70':'8', \
'71':'10', '72':'10', '73':'10', '74':'10', '75':'4', '76':'4', '77':'4', '78':'4', '79':'6', '80':'6', \
'81':'4', '82':'5', '83':'4', '84':'4', '85':'4', '86':'4', '87':'6', '88':'6', '89':'4', '90':'4', \
'91':'4', '92':'4', '93':'4', '94':'4', '95':'4', '96':'4', '97':'6', '98':'6', '99':'4', '100':'4', \
'101':'4', '102':'4', '103':'4', '104':'4', '105':'4', '106':'4', '107':'6', '108':'6', '109':'6', \
'110':'6', '111':'4', '112':'4', '113':'4', '114':'4', '115':'4', '116':'4', '117':'4', '118':'4', '119':'5', '120':'5', \
'121':'5', '122':'5', '123':'4', '124':'4', '125':'4', '126':'4', '127':'4', '128':'4', '129':'4', '130':'4', \
'131':'4', '132':'4', '133':'4', '134':'4', '135':'4', '136':'4', '137':'4', '138':'4', '139':'6', '140':'6', \
'141':'6', '142':'6', '143':'12', '144':'12', '145':'12', '146':'13', '147':'12', '148':'13', '149':'12', '150':'12', \
'151':'12', '152':'12', '153':'12', '154':'12', '155':'13', '156':'12', '157':'12', '158':'12', '159':'12', '160':'13', \
'161':'13', '162':'12', '163':'12', '164':'12', '165':'12', '166':'13', '167':'13', '168':'12', '169':'12', '170':'12', \
'171':'12', '172':'12', '173':'12', '174':'12', '175':'12', '176':'12', '177':'12', '178':'12', '179':'12', '180':'12', \
'181':'12', '182':'12', '183':'12', '184':'12', '185':'12', '186':'12', '187':'12', '188':'12', '189':'12', '190':'12', \
'191':'12', '192':'12', '193':'12', '194':'12', '195':'1', '196':'2', '197':'3', '198':'1', '199':'3', '200':'1', \
'201':'1', '202':'2', '203':'2', '204':'3', '205':'1', '206':'3', '207':'1', '208':'1', '209':'2', '210':'2', \
'211':'3', '212':'1', '213':'1', '214':'3', '215':'1', '216':'2', '217':'3', '218':'1', '219':'2', '220':'3', \
'221':'1', '222':'1', '223':'1', '224':'1', '225':'2', '226':'2', '227':'2', '228':'2', '229':'3', '230':'3'}

	### get input parameters
	chk_error = 0
	fname = cif.split('/')[-1][:-4]
	f = open(cif, 'r')
	lines = f.readlines()
	f.close()
	sp_group = '1'
	i = 0
	while i < len(lines) :
		line = lines[i]
		if "_cell_length_a " in line :
			a = line.strip().split()[1]
			if '(' in a :
				a = a[:a.index('(')]
			if a.replace('.','').isdigit() :
				a = float(a)
			else : chk_error = 1
			i+=1; continue
		if "_cell_length_b " in line :
			b = line.strip().split()[1]
			if '(' in line :
				b = b[:b.index('(')]
			if b.replace('.','').isdigit() :
				b = float(b)
			else : chk_error = 1
			i+=1; continue
		if "_cell_length_c " in line :
			c = line.strip().split()[1]
			if '(' in line :
				c = c[:c.index('(')]
			if c.replace('.','').isdigit() :
				c = float(c)
			else : chk_error = 1
			i+=1; continue
		if "_cell_angle_alpha " in line :
			alpha = line.strip().split()[1]
			if '(' in line :
				alpha = alpha[:alpha.index('(')]
			if alpha.replace('.','').isdigit() :
				alpha = float(alpha)
			else : chk_error = 1
			i+=1; continue
		if "_cell_angle_beta " in line :
			beta = line.strip().split()[1]
			if '(' in line :
				beta = beta[:beta.index('(')]
			if beta.replace('.','').isdigit() :
				beta = float(beta)
			else : chk_error = 1
			i+=1; continue
		if "_cell_angle_gamma " in line :
			gamma = line.strip().split()[1]
			if '(' in line :
				gamma = gamma[:gamma.index('(')]
			if gamma.replace('.','').isdigit() :
				gamma = float(gamma)
			else : chk_error = 1
			i+=1; continue
		if "_symmetry_Int_Tables_number" in line or "_space_group_IT_number" in line :
			sp_group = line.strip().split()[1]
			i+=1; continue
		if "_symmetry_equiv_pos_as_xyz" in line or "_space_group_symop_operation_xyz" in line :
			sym = []
			i+=1; tmp = lines[i].strip()
			while tmp.split()[0].isdigit() :
				tmp = tmp[tmp.index("'"):]
				tmp = tmp.replace("'","").replace("/",".0/")
				sym[len(sym):] = [tmp]
				i+=1; tmp = lines[i]
			continue
		if "_atom_site_occupancy" in line :
			i+=1
			if lines[i].strip()[0] == '_' and not ("_atom_site_attached_hydrogens" in lines[i].strip()) :
				continue

			atoms = []	# Atom postitions
			if "_atom_site_attached_hydrogens" in lines[i] :
				i+=1
			while i < len(lines) :
				tmp = lines[i].split()
				if not(len(tmp) == 7 or len(tmp) == 9 or len(tmp) == 10) :
					break
				if len(tmp) >= 9 :
					for j in [4,5,6,8] :
						if '(' in tmp[j] :
							tmp[j] = float(tmp[j][:tmp[j].index('(')])
						else :
							tmp[j] = float(tmp[j])
						if tmp[j] < 0 :
							tmp[j] = tmp[j]+1
					tmp = [tmp[0]]+tmp[4:7]+[tmp[8],tmp[0],tmp[1]]
					if tmp[0][1].isdigit() :
						tmp[0] = tmp[0][0]
					else :
						tmp[0] = tmp[0][0:2]
				elif len(tmp) == 7 :
					for j in [3,4,5,6] :
						if '(' in tmp[j] :
							tmp[j] = float(tmp[j][:tmp[j].index('(')])
						else :
							tmp[j] = float(tmp[j])
						if tmp[j] < 0 :
							tmp[j] = tmp[j]+1
					tmp = [tmp[0]]+tmp[3:7]+[tmp[1],tmp[0]+'0']
				### Partial occupation process (only for vacancy and multivalency)
				if tmp[4] != 1.0 :
					warn = open(fname+'_vasp2cif_warning!', 'a+')
					warn.write('Check occupancy on '+tmp[5]+'\n')
					chk_occ = 1
					for j in range(len(atoms)) :
						if [tmp[1], tmp[2], tmp[3]] == [atoms[j][1], atoms[j][2], atoms[j][3]] :
							if tmp[0] != atoms[j][0] :
								warn.write('Impossible to decide atom!\n')
								error = open(fname+'_Error','a+')	# ERROR
								error.write('Atom occupancy error on '+cifs[i][:-4])	# ERROR
								chk_error = 1
								break
							else :
								warn.write('Multivalence on '+atoms[j][0]+' - specific magnetic structure cannot be resolved.\n')
								atoms[j][4] = atoms[j][4]+tmp[4]
								chk_occ = 0
								break
					warn.close()
					if chk_occ :
						atoms[len(atoms):]=[tmp]
				### Partail occupation process end
				else :						
					atoms[len(atoms):] = [tmp]
				i+=1
			continue

		if "_atom_site_type_symbol" in line :	# for CoRE MOF DB
			i+=1
			if lines[i].strip()[0] == '_' :
				continue

			atoms = []	# Atom postitions
			while i < len(lines) :
				tmp = lines[i].split()
				if not (len(tmp) == 8) :
					break

				for j in [2,3,4] :
					if '(' in tmp[j] :
						tmp[j] = float(tmp[j][:tmp[j].index('(')])
					else :
						tmp[j] = float(tmp[j])
					if tmp[j] < 0 :
						tmp[j] = tmp[j]+1
				tmp = [tmp[7]]+tmp[2:5]+[float(tmp[6]),tmp[0],tmp[7]+'0']
				### Partial occupation process (only for vacancy and multivalency)
				if tmp[4] != 1.0 :
					warn = open(fname+'_vasp2cif_warning!', 'a+')
					warn.write('Check occupancy on '+tmp[5]+'\n')
					chk_occ = 1
					for j in range(len(atoms)) :
						if [tmp[1], tmp[2], tmp[3]] == [atoms[j][1], atoms[j][2], atoms[j][3]] :
							if tmp[0] != atoms[j][0] :
								warn.write('Impossible to decide atom!\n')
								error = open(fname+'_Error','a+')	# ERROR
								error.write('Atom occupancy error on '+cifs[i][:-4])	# ERROR
								chk_error = 1
								break
							else :
								warn.write('Multivalence on '+atoms[j][0]+' - specific magnetic structure cannot be resolved.\n')
								atoms[j][4] = atoms[j][4]+tmp[4]
								chk_occ = 0
								break
					warn.close()
					if chk_occ :
						atoms[len(atoms):]=[tmp]
				### Partail occupation process end
				else :						
					atoms[len(atoms):] = [tmp]
				i+=1
			continue
		i+=1

	# Double check for partial occupancy.
	for j in range(len(atoms)) :
		if atoms[j][4] != 1.0 :
			warn=open(fname+'_Warnings!', 'a+')
			warn.write('Partial occupancy on '+atoms[j][0]+' atom!')
			warn.close()
	if chk_error == 1 : print("ERROR")	

	table = k_table[sp_group]
	# Calculate lattice vectors
	alpha = math.radians(alpha)
	beta = math.radians(beta)
	gamma = math.radians(gamma)
	if alpha == math.radians(90) and beta ==  alpha and gamma == alpha :
		if table == '8' or table == '9' or table == '10' :
			if a < b and b < c :
				aa = [a,0,0]
				bb = [0,b,0]
				cc = [0,0,c]
				nrepl = 0
			else :
				l_list = [[0,a],[1,b],[2,c]]
				l_list = sorted(l_list, key=itemgetter(1))
				nrepl = 0
				if l_list[0][0] != 0 :
					nrepl = nrepl + 1
				if l_list[1][0] != 1 :
					nrepl = nrepl + 1
				if l_list[2][0] != 2 :
					nrepl = nrepl + 1
				aa = [l_list[0][1],0,0]
				bb = [0,l_list[1][1],0]
				cc = [0,0,l_list[2][1]]
		else :		
			aa = [a,0,0]
			bb = [0,b,0]
			cc = [0,0,c]
	elif table == '15' or table == '16' :
		tmp = a
		a = b
		b = tmp
		tmp = alpha
		alpha = math.radians(180)-beta
		beta = tmp
		aa = [a,0,0]
		bb = [0,b,0]
		cc = [0,c*math.cos(alpha),c*math.sin(alpha)]
	else :
		cg2 = (math.cos(gamma)-math.cos(alpha)*math.cos(beta))/(math.sin(alpha)*math.sin(beta))
		sg2 = (1-cg2**2)**0.5
		cc = [0, 0, c]
		bb = [0, b*math.sin(alpha), b*math.cos(alpha)]
		aa = [a*math.sin(beta)*sg2, a*math.sin(beta)*cg2, a*math.cos(beta)]
		#aa = [0, (a*b*math.cos(gamma)-a*b*math.cos(alpha)*math.cos(beta))/(a*math.sin(alpha)), a*math.cos(beta)]
		aa[0] = math.sqrt(a*a-aa[1]*aa[1]-aa[2]*aa[2])
	### Primitive cell conversion
	if primitive == True :
		if table == '2' or table == '8' or table == '9' :
			[aa,bb,cc] = cvt_cell.lat_fc(aa,bb,cc)
		if table == '3' or table == '5' or table == '6' or table == '10' :
	                [aa,bb,cc] = cvt_cell.lat_bc(aa,bb,cc)
		if table == '11' :
			[aa,bb,cc] = cvt_cell.lat_cc1(aa,bb,cc)
		if table == '13' and gamma == math.radians(120) :
			[aa,bb,cc] = cvt_cell.lat_rh(aa,bb,cc)
#		if table == '16' :
#			[aa,bb,cc] = cvt_cell.lat_cc2(aa,bb,cc)

	### Symmetry reproduction
	atom_z = []		# Final atom names
	atom_cnt = []	# Final atom counts
	final = []
	k = 0
	chem = atoms[0][0]
	atom_z.append(chem)
	atom_cnt.append(0)
	for index in range(len(atoms)) :
		if chem != atoms[index][0] :
			chem = atoms[index][0]
			atom_z.append(chem)
			atom_cnt.append(0)
			k = k+1
		# final[i] = [x,y,z,atom_type,atom_oxi]
		final[len(final):] = [cvt_cell.setBasis(atoms[index][1],atoms[index][2],atoms[index][3])]
		final[len(final)-1][3:] = [atoms[index][5], atoms[index][6]]
		atom_cnt[k] = atom_cnt[k]+1
		for j in range(len(sym)) :
			[x, y, z] = [atoms[index][1], atoms[index][2], atoms[index][3]]
			xx = eval(sym[j].split(', ')[0])
			yy = eval(sym[j].split(', ')[1])
			zz = eval(sym[j].split(', ')[2])
			[xx,yy,zz] = cvt_cell.setBasis(xx,yy,zz)
			check = 0
			for l in range(len(final)) :
				if abs(xx-final[l][0]) < 0.001 and abs(yy-final[l][1]) < 0.001 and abs(zz-final[l][2]) < 0.001 :
					check = 1
			if check == 0 :
				final[len(final):] = [[xx, yy, zz, atoms[index][5], atoms[index][6]]]
				atom_cnt[k] = atom_cnt[k]+1
	if table == '8' or table =='9' or table == '10' :
		if nrepl == 2 :
			for j in range(len(final)) :
				tt = [final[j][0],final[j][1],final[j][2]]
				[final[j][0],final[j][1],final[j][2]] = [1-tt[l_list[0][0]],tt[l_list[1][0]],tt[l_list[2][0]]]
		if nrepl == 3 :
			for j in range(len(final)) :
				tt = [final[j][0],final[j][1],final[j][2]]
				[final[j][0],final[j][1],final[j][2]] = [tt[l_list[0][0]],tt[l_list[1][0]],tt[l_list[2][0]]]
	if table == '15' or table == '16' :
		for j in range(len(final)) :
			[t1,t2,t3] = [final[j][0],final[j][1],final[j][2]]
			[final[j][0],final[j][1],final[j][2]] = [t2,1.0-t1,t3]

	# Save POSCAR object
	name = fname
	sym = sp_group
	mag = 1.0
	lattice = [aa,bb,cc]
	element = atom_z
	natom = atom_cnt
	SD = True
	coord = 'Direct'
	xyz = []; xyz_sd = []; xyz_tag = []
	for i in range(len(final)) :
		xyz.append(final[i][0:3])
		xyz_sd.append(['T','T','T'])
		xyz_tag.append('!'+final[i][3])
	return POSCAR(name,mag,lattice,element,natom,SD,coord,xyz,xyz_sd,xyz_tag,sym)

### Main ###
if __name__ == '__main__':
	cifs = sys.argv[1]
	poscar_out = cif2vasp(cifs,primitive=False)
	poscar_out.merge_element()
	fname = poscar_out.name
	out = open('POSCAR_'+fname, 'w')
	out.write(poscar_out.make_POSCAR())
	out.close()
