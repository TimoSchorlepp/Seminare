import numpy as np
import matplotlib.pyplot as plt
import sys

np.random.seed(42)

def init_config(m,init_type='zero',bc='nonperiodic'):
	ret = np.zeros((m,m))
	if init_type == 'zero':
		pass
	elif init_type == 'full':
		if bc == 'nonperiodic':
			ret[::2,::2] = np.ones((m/2,m/2))
			ret[1::2,1::2] = np.ones((m/2,m/2))
		elif bc == 'periodic':
			ret[2::2,2::2] = np.ones((m/2-1,m/2-1))
			ret[1:-2:2,1:-2:2] = np.ones((m/2-1,m/2-1))
		else:
			print 'Error: Unknown boundary condition requested'
	else:
		print 'Error: Unknown initial configuration requested'
	return ret

def get_random_pos(m):
	return np.random.randint(0,m),np.random.randint(0,m)

def is_valid_flip(i,j,m,occ_pos,bc='nonperiodic'):
	ret = True
	if occ_pos[i,j] == 1:
		pass
	else:
		if bc == 'nonperiodic':
			n1 = (i,j) if j == m-1 else (i,j+1)
			n2 = (i,j) if j == 0 else (i,j-1)
			n3 = (i,j) if i == m-1 else (i+1,j)
			n4 = (i,j) if i == 0 else (i-1,j)
		elif bc == 'periodic':
			n1 = (i,(j+1)%m)
			n2 = (i,(j-1)%m)
			n3 = ((i+1)%m,j)
			n4 = ((i-1)%m,j)
		else:
			print 'Error: Unknown boundary condition requested'
		if occ_pos[n1]+occ_pos[n2]+occ_pos[n3]+occ_pos[n4]>0:
			ret = False
	return ret

def flip(i,j,occ_pos):
	ret = occ_pos
	ret[i,j] = 1 - ret[i,j]
	return ret

def get_rel_occupancy(occ_pos):
	return np.mean(occ_pos)

def MCMC(m,init_type,bc,burn_in,N):
	occ_pos = init_config(m,init_type,bc)
	
	for k in range(burn_in):
		i,j = get_random_pos(m)
		if is_valid_flip(i,j,m,occ_pos,bc):
			occ_pos = flip(i,j,occ_pos)
	
	rel_occupancy_list = np.zeros(N)
	
	for k in range(N):
		sys.stdout.write('Step {} of {}\r'.format(k,N))
		sys.stdout.flush()
		i,j = get_random_pos(m)
		if is_valid_flip(i,j,m,occ_pos,bc):
			occ_pos = flip(i,j,occ_pos)
		rel_occupancy_list[k] = get_rel_occupancy(occ_pos)
	print ' '
	return rel_occupancy_list

def get_data(M,N,burn_in):
	for m in M:
		res = MCMC(m,'zero','nonperiodic',burn_in,N)
		np.save('data/{}_zero_np.npy'.format(m),res)
		res = MCMC(m,'full','nonperiodic',burn_in,N)
		np.save('data/{}_full_np.npy'.format(m),res)
		res = MCMC(m,'zero','periodic',burn_in,N)
		np.save('data/{}_zero_p.npy'.format(m),res)
		res = MCMC(m,'full','periodic',burn_in,N)
		np.save('data/{}_full_p.npy'.format(m),res)

def plot_data(M,N):
	bc = 'np'
	plt.figure()
	N_arr = np.linspace(1,N,N)
	for m in M:
		#~ print m
		res1 = np.load('data/{}_zero_{}.npy'.format(m,bc))
		res2 = np.load('data/{}_full_{}.npy'.format(m,bc))
		mn1 = np.cumsum(res1)/N_arr
		mn2 = np.cumsum(res2)/N_arr
		#~ start = int(9e6)
		#~ print np.mean(mn1[start:]), '+-', np.std(mn1[start:])
		#~ print np.mean(mn2[start:]), '+-', np.std(mn2[start:])
		p = plt.plot(N_arr[::4],mn1[::4],label=r'$m = {{{}}}$'.format(m))
		plt.plot(N_arr[::4],mn2[::4],color=p[0].get_color())
	plt.semilogx()
	plt.legend(loc='best')
	plt.xlabel(r'$n$')
	plt.ylabel(r'$M_n(f)(\omega)$')
	plt.show()
	plt.savefig('hardcore.pdf')

def animate(m,N,interval,init_type,bc):
	occ_pos = init_config(m,init_type,bc)
	mn = 0.
	fig,ax = plt.subplots(1,2)
	ax[1].set_title('Mean occupancy fraction')
	
	for k in range(N):
		i,j = get_random_pos(m)
		if is_valid_flip(i,j,m,occ_pos,bc):
			occ_pos = flip(i,j,occ_pos)
		mn += get_rel_occupancy(occ_pos)
		if k % interval == 0:
			p = ax[0].imshow(occ_pos,cmap='plasma')
			cb = plt.colorbar(p,ax=ax[0])
			ax[0].set_title(r'Configuration at $N = {{{}}}$ steps'.format(k))
			ax[1].scatter(k,mn/k,color='red')
			plt.tight_layout()
			plt.pause(0.01)
			cb.remove()
			ax[0].clear()
			
if __name__ == '__main__':
	M = [8,16,32,64,128]
	N = int(1e7)
	burn_in = 0
	
	#~ get_data(M,N,burn_in)
	#~ plot_data(M,N)
	animate(128,100000,500,'zero','periodic')
		
