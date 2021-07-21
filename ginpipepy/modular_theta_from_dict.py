#modular_ttheta_from_dict.py

"""
Created on Sun Jul 12 11:56:45 2020

@author: Maria Trofimova

Based on: Bhavin S Khatri, Austin Burt, Robust Estimation of Recent
Effective Population Size from Number of Independent Origins in Soft
Sweeps, Molecular Biology and Evolution,
Volume 36, Issue 9, September 2019,
Pages 2040–2052, https://doi.org/10.1093/molbev/msz081
"""

"""
Metrics from bins
"""
import numpy as np
import math
from scipy import optimize


class analyzeTrajectory:

    def __init__(self, dict_traj, mut_proportion, num_days_per_bin, initSeq):
        '''Class: calculate point-wise constant effective population size estimate
        :param dict_traj: list of sequence fingerprints in format ('id','MUT_POS_1>MUT_BASE-MUT_POS_2>MUT_BASE-...')  
        :type dict_traj: list
        :param mut_proportion: proportion of mutant bases in entire data set
        :type mut_proportion: float
        :param num_days_per_bin: list of number of days hat each bin in dict_traj is spanning
        :type num_days_per_bin: list[int]
        :param initSeq: initial sequence
        :type initSeq: str
        '''
        self.dict_traj = dict_traj
        self.initSeq = initSeq
        self.mut_proportion = mut_proportion
        self.delta_t = num_days_per_bin

    def _makeOriginsFirstOcc(self, mutSeq, mutantsCount):
        '''
        Find first occurence of mutant - origins from sequences not recorded
        before current bin
        :param mutSeq: dictionary with key - sequence, value - list of length
        t - number of bins, with 1 if sequence present in the bin and 0 otherwise
        :type mutSeq: dict
        :param mutantsCount: list of length t - number of bins, with
        number of mutant sequences in bin
        :type mutantsCount: list
        :returns: list of length t - number of unique first occuring sequences in bin
        :rtype: list
        '''
        origins = []
        firstOcc = np.zeros(len(mutSeq))
        for i,(key,val) in enumerate(mutSeq.items()):
            # Index of first occurence
            ind = np.where(val==1)[0][0]
            firstOcc[i] = ind

        for i in range(len(mutantsCount)):
            # First occurences in current bin
            c = len(np.where(firstOcc==(i))[0])
            origins.append(c)

        return origins

    def _findPosition(self, index, occList):
        # List of indices of lists where the index was found
        bin_indices = []
        for i, pos in enumerate(occList):
            val_list = np.where(pos==index)[0]
            if val_list.size!=0:
                bin_indices.append(i)
        return bin_indices

    def _makeOrigins(self, mutSeq, mutantsCount):
        '''
        Find all occurences of mutant - count unique mutant types in current bin
        :param mutSeq: dictionary with key - sequence, value - list of length
        t - number of bins, with 1 if sequence present in the bin and 0 otherwise
        :param mutantsCount: list of lengtth t - number of bins, with
        number of mutant sequences in bin
        :type mutSeq: dict
        :returns: list of length t - number of unique sequences in bin
        :rtype: list
        '''
        origins = []
        occurence = [[] for i in range(len(mutSeq))]
        for i,(key,val) in enumerate(mutSeq.items()):
            ind = np.where(val==1)[0]
            occurence[i] = ind

        for i in range(len(mutantsCount)):
            ind_list = self._findPosition(i,occurence)
            origins.append(len(ind_list))

        return origins

    def _fmle(self,x,nu,ns):
        '''
        Function of the estimate for root finding
        :param nu: number of origins
        :type nu: int
        :param ns: number of mutant sequences
        :type ns: int
        :param x: param to optimise
        :type x: float
        :returns: population size estimate
        :rtype: float
        '''
        a = float(x*np.log(1+ns/x))
        b = np.math.factorial(nu)
        #c = np.exp(-(x*np.log(1+ns/x)))
        d = nu*math.log(a)
        dd = math.log(b)
        e = x*math.log(1+ns/x)
        return -(d-dd-e)

    def _fmle_short(self,x,nu,ns):
        '''
        Function of the estimate for root finding
        :param nu: number of origins
        :type nu: int
        :param ns: number of mutant sequences
        :type ns: int
        :param x: param to optimise
        :type x: float
        :returns: population size estimate
        :rtype: float
        '''
        return (nu-(x*np.log(1+ns/x)))**2

    def _optimize(self, nu, ns):
        """
        Optimization routines
        :param nu: number of origins
        :type nu: int
        :param ns: number of mutant sequences
        :type ns: int
        :returns: optimal population size estimate
        :rtype: float
        """
        # Constraint
        con1 = {'type': 'ineq', 'fun': lambda x: x}
        # If no data available -- set estimate to zero
        if nu==0 or ns==0:
            return 0
        # Else optimize with trust-constr option to keep evaluation in feasible region
        else:
            x0 = 1
            #sol = optimize.minimize(self._fmle, x0=x0, method='trust-constr', args=(nu,ns), constraints=con1)
            # use least squares without MLE
            sol = optimize.minimize(self._fmle_short, x0=x0, method='trust-constr', args=(nu,ns), constraints=con1)
            return sol.x[0]

    def analyzeBinsMLE(self):
        """
        Estimates joint parameter 2Nmu from data
        Count mutant sequences per bin and creates sequence dictionary with
        bin indices where the mutant was found.
        :param:
        :returns: thetas list - parameter estimates; list variance - variance estimates
        from maximum likelihood estimate, variance_size - 1/size_of_bin, num_mut - number of
        mutants per bin, num_seqs - number of sequences per bin, origins - number
        of origins per bin
        :rtype: list
        """
        # Number of (mutated) sequences
        num_mut = []
        num_seqs = []
        # Dictionary of mutants
        # key: unique sequence
        # value: list of length t with 1 if sequence presentt in bin t and
        # 0 otherwise
        mutSeqDict = dict()

        # Length of sequence evolution
        traj_length = len(self.dict_traj)

        # Filling the mutants dict and presence array
        # Variance from bin size - 1/bin_size
        variance_size = []
        for t, seqSet in enumerate(self.dict_traj):
            mut_count = 0
            # If bin is not empty
            if len(seqSet)!=0:
                #print("Sample size: ",len(seqSet))
                num_seqs.append(len(seqSet))
                variance_size.append(1/len(seqSet))
                for i in range(len(seqSet)):
                    # If the sequence has mutant bases
                    if (seqSet[i][1]!=''):
                        mut_count += 1
                        # Fill the mutant dictionary
                        # If the sequence FP was not seen before
                        if not seqSet[i][1] in mutSeqDict:
                            # Zero filled list of length n_bins
                            times = np.zeros(traj_length)
                            # 1 for current bin
                            times[t] = 1
                            mutSeqDict[seqSet[i][1]] = times
                        # If the sequence FP was seen before
                        else:
                            times = mutSeqDict[seqSet[i][1]]
                            # If already seen in current bin - do nothing, else -
                            # replace the current bin with 1
                            if times[t]==1:
                                pass
                            elif times[t]==0:
                                times[t] = 1
                                mutSeqDict[seqSet[i][1]] = times
                num_mut.append(mut_count)
            # If bin empty
            else:
                num_seqs.append(0)
                variance_size.append(1)
                num_mut.append(mut_count)


        origins = self._makeOrigins(mutSeqDict, num_mut)
        #Estimate effective population size from data
        thetas = []
        #Variance
        variance = []

        for i in range(len(origins)):
            # Get the MLE
            #sol = self._optimize(origins[i], num_mut[i])
            # Get direct fit original
            sol = self._optimize(origins[i]/(1+np.log(math.sqrt(self.delta_t[i]))), num_mut[i]/(1+np.log(math.sqrt(self.delta_t[i]))))
            thetas.append(sol)

        for i in range(len(thetas)):
            # Variance of the estimate
            if thetas[i]!=0:
                var = thetas[i]*math.log(1+num_mut[i]/thetas[i])
                variance.append(var)
            else:
                variance.append(0)

        return thetas, variance, variance_size, num_seqs, num_mut, origins
