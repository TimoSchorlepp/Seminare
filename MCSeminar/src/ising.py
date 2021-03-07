import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)
beta_c = 0.5*np.log(1+np.sqrt(2)) # \approx 0.44069

def init_spins(m,init_type='unmagnetized'):
	ret = np.ones((m,m))
	if init_type == 'magnetized':
		pass
	elif init_type == 'unmagnetized':
		ret = 2*(np.random.randint(0,2,size=(m,m)) - 0.5)
	else:
		print 'Error: Unknown initial configuration requested'
	return ret

def get_initial_energy(spins,m):
	ener = 0
	for i in range(m):
		for j in range(m):
			sp1 = spins[(i+1)%m,j]
			sp2 = spins[i,(j+1)%m]
			ener += - spins[i,j] * (sp1+sp2)
	return ener

def update_energy(ener,de):
	return ener+delta

def get_random_pos(m):
	return np.random.randint(0,m),np.random.randint(0,m)

def get_rel_magnetization(spins,abs_val=False):
	if abs_val:
		return abs(np.mean(spins))
	else:
		return np.mean(spins)

def get_energy_difference(spins,m,i,j):
	sp1 = spins[(i+1)%m,j]
	sp2 = spins[(i-1)%m,j]
	sp3 = spins[i,(j+1)%m]
	sp4 = spins[i,(j-1)%m]
	return 2 * spins[i,j] * (sp1+sp2+sp3+sp4)

def get_acceptance_rate(beta,ener,de,sampler='metropolis'):
	ret = 0
	if sampler == 'metropolis':
		ret = np.exp(-beta*de)
	elif sampler == 'gibbs':
		ret = np.exp(-beta*(ener+de))/(np.exp(-beta*ener)+np.exp(-beta*(ener+de)))
	else:
		print 'Error: Unknown sampler requested'
	return ret

def flip_spin(spins,i,j):
	ret = spins
	ret[i,j] = -ret[i,j]
	return ret

def MCMC(m,beta,init_type,burn_in,N,sampler):
	spins = init_spins(m,init_type)
	ener = get_initial_energy(spins,m)
	
	for k in range(burn_in):
		i,j = get_random_pos(m)
		de = get_energy_difference(spins,m,i,j)
		alpha = get_acceptance_rate(beta,ener,de,sampler)
		if np.random.uniform() <= alpha:
			spins = flip_spin(spins,i,j)
			ener = ener+de
		if k % (burn_in/100) == 0:
			print k*100./burn_in
	
	magnetization_list = np.zeros(N)
	
	for k in range(N):
		i,j = get_random_pos(m)
		de = get_energy_difference(spins,m,i,j)
		alpha = get_acceptance_rate(beta,ener,de,sampler)
		if np.random.uniform() <= alpha:
			spins = flip_spin(spins,i,j)
			ener = ener+de
		magnetization_list[k] = get_rel_magnetization(spins)
		if k % (N/100) == 0:
			print k*100./N
	
	return magnetization_list
	

def animate(m,beta,N,interval,init_type,sampler):
	spins = init_spins(m,init_type)
	ener = get_initial_energy(spins,m)
	magnetization_cumulative = get_rel_magnetization(spins)
	ener_cumulative = ener
	fig,ax = plt.subplots(1,3)
	ax[1].set_title(r'$<M>$')
	ax[2].set_title(r'$<E>$')
	
	for k in range(N):
		i,j = get_random_pos(m)
		de = get_energy_difference(spins,m,i,j)
		alpha = get_acceptance_rate(beta,ener,de,sampler)
		if np.random.uniform() <= alpha:
			spins = flip_spin(spins,i,j)
			ener = ener+de
		magnetization_cumulative += get_rel_magnetization(spins,abs_val=False)
		ener_cumulative += ener
		if k % interval == 0:
			p = ax[0].imshow(spins,cmap='plasma')
			cb = plt.colorbar(p,ax=ax[0])
			ax[0].set_title(r'Configuration at $N = {{{}}}$ steps'.format(k))
			ax[1].scatter(k,magnetization_cumulative/k,color='blue')
			ax[2].scatter(k,ener_cumulative/k/m**2,color='red')
			plt.tight_layout()
			plt.pause(0.01)
			cb.remove()
			ax[0].clear()

def get_data_magnetization_diff_beta(M,beta,N,burn_in):
	for bet in beta:
		for m in M:
			mag = MCMC(m,bet,'unmagnetized',burn_in,N,'metropolis')
			np.save('data/{}_{}_mag_un.npy'.format(bet,m),mag)
			mag = MCMC(m,bet,'magnetized',burn_in,N,'metropolis')
			np.save('data/{}_{}_mag_full.npy'.format(bet,m),mag)

def get_data_abs_magnetization_diff_beta_diff_res(M,beta,N,burn_in):
	for bet in beta:
		for m in M:
			mag = MCMC(m,bet,'unmagnetized',burn_in,N,'metropolis')
			np.save('data/abs_{}_{}_mag.npy'.format(bet,m),np.abs(mag))


def plot_data_magnetization_diff_beta(M,beta,N):
	plt.figure()
	N_arr = np.linspace(1,N,N)
	for bet in beta:
		print bet
		for m in M:
			res1 = np.load('data/{}_{}_mag_un.npy'.format(bet,m))
			res2 = np.load('data/{}_{}_mag_full.npy'.format(bet,m))
			mn1 = np.cumsum(res1)/N_arr
			mn2 = np.cumsum(res2)/N_arr
			#~ start = int(9e6)
			#~ print np.mean(mn1[start:]), '+-', np.std(mn1[start:])
			#~ print np.mean(mn2[start:]), '+-', np.std(mn2[start:])
			p = plt.plot(N_arr[::20],mn1[::20],label=r'$\beta = {{{}}}$'.format(round(bet,4)))
			plt.plot(N_arr[::20],mn2[::20],color=p[0].get_color())
	plt.semilogx()
	plt.legend(loc='best')
	plt.xlabel(r'$n$')
	plt.ylabel(r'$M_n(f)(\omega)$')
	plt.show()
	#~ plt.savefig('ising_magnetization_250_init.pdf')

def plot_data_abs_magnetization_diff_beta_diff_res(M,beta,N):
	plt.figure()
	N_arr = np.linspace(1,N,N)
	clist = ['blue','orange','green']
	counter = 0
	# ~ for m in M:
		# ~ plt.scatter([],[],color=clist[counter],label=r'$m = {{{}}}$'.format(m))
		# ~ for bet in beta:
			#~ print m
			# ~ res1 = np.load('data/abs_{}_{}_mag.npy'.format(bet,m))
			# ~ mn1 = np.cumsum(res1)/N_arr
			#~ start = int(9e6)
			#~ print np.mean(mn1[start:]), '+-', np.std(mn1[start:])
			#~ print np.mean(mn2[start:]), '+-', np.std(mn2[start:])
			# ~ p = plt.scatter(bet,mn1[-1],color=clist[counter])
		# ~ counter += 1
	# plot exact result for m \to \infty:
	n = 500
	beta_plt = np.linspace(0.1,0.9,n)
	Ef = np.zeros(n)
	for i in range(n):
		if beta_plt[i]>beta_c:
			Ef[i] = (1-np.sinh(2*beta_plt[i])**(-4))**(1./8.)
	plt.plot(beta_plt,Ef,c='red')
	# ~ plt.legend(loc='best')
	plt.xlabel(r'$\beta$')
	# ~ plt.ylabel(r'$M_n(f)(\omega)$')
	plt.ylabel(r'$E_{\mu^\beta}(f)$')
	plt.show()
	#~ plt.savefig('ising_magnetization_250_init.pdf')

if __name__ == '__main__':
	# ~ M = [128,256,1024]
	# ~ beta = [0.2,0.25,0.3,0.35,0.4,beta_c,0.5,0.55,0.6,0.7,0.8]
	# ~ N = int(1e6)
	# ~ burn_in = 10*N
	# ~ get_data_abs_magnetization_diff_beta_diff_res(M,beta,N,burn_in)
	# ~ plot_data_abs_magnetization_diff_beta_diff_res(M,beta,N)
	
	# ~ M = [256]
	# ~ beta = [0.2,0.4,beta_c,0.6]
	# ~ N = int(1e8)
	# ~ burn_in = 0
	# ~ get_data_magnetization_diff_beta(M,beta,N,burn_in)
	# ~ plot_data_magnetization_diff_beta(M,beta,N)
	
	animate(64,100*beta_c,10000000,2000,'unmagnetized','metropolis')
