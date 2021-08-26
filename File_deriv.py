import numpy as np
import CAMB_base as cb
import camb
import os 

#fonction de bruit
def noise(s,l,theta_FWHM):
    return s**2*np.exp(l*(l+1)*theta_FWHM**2/(8*np.log(2)))


#nous renvoie la fonction en escalier
def construction(bins):
    taille = np.shape(bins)[0]+6
    z = np.zeros(taille)
    z[0], z[1],z[2],z[3]= 1,3.8,4.3,6
    z[taille-2],z[taille-1]=30.25,35
    for i in range(4,taille-2):
        z[i]=6+(i-3)*0.25
    xe = np.zeros(np.shape(z)[0])
    xe[0],xe[1],xe[2],xe[3]= 1.08,1.08,1,1
    xe[taille-2],xe[taille-1]=0,0
    for i in range(4,taille-2):
        xe[i] = bins[i-4]
    return z,xe
    
    #z = np.array([1,2.9999,3,3.9999,4,4.9999,5,5.9999,6,6.9999,7,7.9999,8,8.9999,9,9.9999,10,15])
    xe = np.zeros(np.shape(z)[0])
    xe[0],xe[1],xe[2],xe[3]= 1.08,1.08,1,1
    taille = np.shape(z)[0]
    xe[taille-1],xe[taille-2]=0,0
    for i in range(2,taille//2-1):
        xe[2*i],xe[2*i+1]= bins[i-2],bins[i-2]
    return z,xe



#calcule la dérivée par rapport à un bin
def get_deriv_bin(parameter,liste,bins,amplitude):
    cmb_params = {"H0": 67.4,"As": 1e-9,"ombh2": 0.02237,"omch2": 0.1200,"ns": 0.9649 }
    
    if parameter not in liste:
        if parameter == 'As':
            stepsize = 10**(-11)
        if parameter != 'As':
            stepsize = 0.01
        
        cmb_params_l = cmb_params.copy()
        cmb_params_r = cmb_params.copy()
        model = camb.TanhReionization()
        z,xe = construction(bins)
        model.set_xe(z,xe,smooth=0)
        cmb_params_l[parameter] = cmb_params[parameter]-stepsize
        cmb_params_r[parameter] = cmb_params[parameter]+stepsize
        pars_left = camb.set_params(**cmb_params_l,Reion=model)
        pars_right = camb.set_params(**cmb_params_r,Reion=model)
        
    if parameter in liste :
        stepsize=0.001
        bins_l = bins.copy()
        bins_r = bins.copy()
        bins_l[amplitude[parameter]-1] =   bins_l[amplitude[parameter]-1]-stepsize
        bins_r[amplitude[parameter]-1] = bins_r[amplitude[parameter]-1]+ stepsize
        model_l = camb.TanhReionization()
        model_r = camb.TanhReionization()
        z_l,xe_l = construction(bins_l)
        z_r,xe_r = construction(bins_r)
        model_l.set_xe(z_l,xe_l,smooth=0)
        model_r.set_xe(z_r,xe_r,smooth=0)
        pars_left = camb.set_params(**cmb_params,Reion=model_l)
        pars_right = camb.set_params(**cmb_params,Reion=model_r)
        #permet de modifier la précision de CAMB
        pars_left.set_accuracy(AccuracyBoost=7.0, lSampleBoost=7.0, lAccuracyBoost=7.0)
        pars_right.set_accuracy(AccuracyBoost=7.0, lSampleBoost=7.0, lAccuracyBoost=7.0)
        

    result_left = camb.get_results(pars_left)
    result_right = camb.get_results(pars_right)
    powers_left =result_left.get_cmb_power_spectra(pars_left,raw_cl=True,lmax=500)
    powers_right =result_right.get_cmb_power_spectra(pars_right,raw_cl=True,lmax=500)   
    cl_left = powers_left['total']
    cl_right = powers_right['total']
    cl_tt_left = cl_left[:,0]
    cl_tt_right = cl_right[:,0]
    cl_te_left = cl_left[:,3]
    cl_te_right = cl_right[:,3]
    cl_ee_left = cl_left[:,1]
    cl_ee_right = cl_right[:,1]
    cl_bb_left = cl_left[:,2]
    cl_bb_right = cl_right[:,2]
       
    dCltt_dh = (cl_tt_right - cl_tt_left) / (2 * stepsize)
    dClte_dh = (cl_te_right - cl_te_left) / (2 * stepsize)
    dClee_dh = (cl_ee_right - cl_ee_left) / (2 * stepsize)
    dClbb_dh = (cl_bb_right - cl_bb_left) / (2 * stepsize)
    return dCltt_dh,dClte_dh,dClee_dh,dClbb_dh



    
sT = 2.2 * (np.pi/60./180.)*10**(-6)*0
sP = sT * np.sqrt(2.)*0
theta_FWHM = 30. * (np.pi/60./180.)*0
f_sky = 1

liste=[]
for i in range(1,97):
    liste.append('a'+str(i))
amplitude={}
for i in range(1,97):
    amplitude["a"+str(i)]=i
    
bins = np.full(96,0.75)

#permet d'obtenir toutes les dérivées par rapport à un bin
arrayID = os.environ["SGE_TASK_ID"]
parameter = 'a'+str(arrayID)
Deriv = get_deriv_bin(parameter,liste,bins,amplitude)[2]
np.savetxt('/pbs/home/e/epangbur/stagem1/deriv_10_'+str(arrayID)+'.txt',Deriv)


#cette partie permet d'obtenir le c_l associé à une fonction test
"""
cmb_params = {"H0": 67.4,"As": 1e-9,"ombh2": 0.02237,"omch2": 0.1200,"ns": 0.9649 }
model = camb.TanhReionization()
z,xe = construction(bins)
model.set_xe(z,xe,smooth=0)   
pars = camb.set_params(**cmb_params,Reion=model)
pars.set_accuracy(AccuracyBoost=7.0, lSampleBoost=7.0, lAccuracyBoost=7.0)
results = camb.get_results(pars)
powers =results.get_cmb_power_spectra(pars,raw_cl=True,lmax=500)
totCl=powers['total'][:,2]
np.savetxt('/pbs/home/e/epangbur/stagem1/c_l10.txt',totCl)
"""