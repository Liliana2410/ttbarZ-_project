import ROOT
from ROOT import TLorentzVector, TH1F
import histos
from collections import Counter
from math import sin, cos, sinh


def jet_reconstruction(jets, bs):

	global dijet_1
	global b_dijet_1
	global dijet_2
	global b_dijet_2
	global reconstructed_W1
	global reconstructed_W2
	#merged category of the first pair of hadronic jets coming from the top.
	global notMerged_1
	global partiallyMerged_1
	global fullyMerged_1
	#merged category of the second pair of hadronic jets coming from the top.
	global notMerged_2
	global partiallyMerged_2
	global fullyMerged_2
	#index of the first b-jet from the hadronic decaying top.
	global index_b1
	#index of the second b-jet from the hadronic decaying top.
	global index_b2
	#index of the two pairs of jets coming from the hadronic tops (pairs are always j1-j2 and j3-j4).
	global index_j1
	global index_j2
	global index_j3
	global index_j4

	#global Delta_PT_b_alpha_b_beta

	best_Err = 999999.9
	MW = 80.379
	Mt = 172.76

	index_j_temp1 = 0
	index_j_temp2 = 0
	index_j_temp3 = 0
	index_j_temp4 = 0
	
	index_b_temp1 = 0
	index_b_temp2 = 0
	
	#loop over the first group of particles (j_1, j_2, b_1).
	for i in range(len(jets)):
		for j in range(i + 1, len(jets)):
				for k in range(len(bs)):
				
					#Loop over the second group of particles (j_3, j_4, b_2).
					for n in range(len(jets)):
						if ((n != i) and (n != j)):
							for m in range(n + 1, len(jets)):
								if ((m != i) and (m != j)):
									for l in range(len(bs)):
										if (l != k):
							
											jtemp1 = jets[i]
											jtemp2 = jets[j]
											jtemp3 = jets[n]
											jtemp4 = jets[m]
											btemp1 = bs[k]
											btemp2 = bs[l]
												
											Err_1 = (abs((jtemp1 + jtemp2 + btemp1).M() - Mt))*100/Mt
											Err_2 = (abs((jtemp3 + jtemp4 + btemp2).M() - Mt))*100/Mt
											Err = Err_1 + Err_2
					
											#Selection criteria.
											if (Err < best_Err):
												best_Err = Err
												index_j_temp1 = i
												index_j_temp2 = j
												index_j_temp3 = n
												index_j_temp4 = m
												index_b_temp1 = k
												index_b_temp2 = l
					
	
	#Initialize the first group of particles.				
	index_j1 = index_j_temp1
	index_j2 = index_j_temp2	
	index_b1 = index_b_temp1
	
	dijet_1 = jets[index_j1] + jets[index_j2]		
	dr_dijet_1 = jets[index_j1].DeltaR(jets[index_j2])
	b_dijet_1 = jets[index_j1] + jets[index_j2] + bs[index_b1]
	dr_b_dijet_1 = dijet_1.DeltaR(bs[index_b1])
	
	#Initialize the second group of particles.
	index_j3 = index_j_temp3
	index_j4 = index_j_temp4	
	index_b2 = index_b_temp2
	
	dijet_2 = jets[index_j3] + jets[index_j4]		
	dr_dijet_2 = jets[index_j3].DeltaR(jets[index_j4])
	b_dijet_2 = jets[index_j3] + jets[index_j4] + bs[index_b2]
	dr_b_dijet_2 = dijet_2.DeltaR(bs[index_b2])
    		
	#Merged category for the first group of particles.
	if (dr_dijet_1 > 0.8):

	 	notMerged_1 = True
	 	partiallyMerged_1 = False
	 	fullyMerged_1 = False


   	else:

		notMerged_1 = False
		
        	if (dr_b_dijet_1 > 1.0):
        
			partiallyMerged_1 = True
			fullyMerged_1 = False
			reconstructed_W1 = dijet_1

		else:

			partiallyMerged_1 = False
			fullyMerged_1 = True
		
	#Merged category for the second group of particles.	
	if (dr_dijet_2 > 0.8):

	 	notMerged_2 = True
	 	partiallyMerged_2 = False
	 	fullyMerged_2 = False


   	else:

		notMerged_2 = False
		
        	if (dr_b_dijet_2 > 1.0):
        
			partiallyMerged_2 = True
			fullyMerged_2 = False
			reconstructed_W2 = dijet_2

		else:

			partiallyMerged_2 = False
			fullyMerged_2 = True

	return jets[index_j1], jets[index_j2], jets[index_j3], jets[index_j4], bs[index_b1], bs[index_b2]


def tau_reconstruction(taus):

	#index of the two taus (pairs of jets tagged as coming from taus), does NOT include the energy from the neutrinos (MET).
	global index_tau1
	global index_tau2

	best_dPt = 99999999.9
	for i in range(len(taus)):
		for j in range(i + 1, len(taus)):
			
			dPt = abs(taus[i].Pt() - taus[j].Pt())
			if (dPt < best_dPt):
			
				best_dPt = dPt
				index_tau1 = i
				index_tau2 = j
					
	return taus[index_tau1], taus[index_tau2]


def cross_cleaning(jets, bs, taus):

	particles = jets + bs
	temp = []
	for i in range(len(taus)):
		
		crossed = False
		for j in range(len(particles)):

			if (taus[i].DeltaR(particles[j]) < 0.3):
				crossed = True

		if (crossed == False):
			
			temp.append(taus[i])

	return temp
	
	
def frac_twobody_pTvspL(particle1, particle2):

	px1 = (particle1.Pt())*(cos(particle1.Phi()))
	py1 = (particle1.Pt())*(sin(particle1.Phi()))
	pz1 = (particle1.Pt())*(sinh(particle1.Eta()))
	
	px2 = (particle2.Pt())*(cos(particle2.Phi()))
	py2 = (particle2.Pt())*(sin(particle2.Phi()))
	pz2 = (particle2.Pt())*(sinh(particle2.Eta()))

	pT = ((px1 + px2)**2 + (py1 + py2)**2)**0.5
	pL = abs(pz1 + pz2)
	
	p = ((px1 + px2)**2 + (py1 + py2)**2 + (pz1 + pz2)**2)**0.5

	frac = pT/p

	return frac


def PT(TLV):
	return TLV.Pt()

def P(TLV):
	return TLV.P()
    
def histos_fill(plot, variable):
	return plot.Fill(variable)


def histos_Draw(plot):
	return plot.Draw('HISTOS')  


def histos_Write(plot):
	return plot.Write()


def histos_Reset(plot):
	return plot.Reset('ICESM') 


entries = ["Zprime_tata_350", "Zprime_tata_1000", "Zprime_tata_1500", "Zprime_tata_3000", "ttbarh", "ttbarZ"]
jobs = [1,1,1,1,1,1]


#Boolean value that determines if plots will be generated or not.
makeplots = False
#Boolean value that determines if the taus should be checked (amount of misidentified taus and tau2 distributions).
checktaus = False
#Boolean value that determines if the gen studies are carried out or not (matching algorithm, etc).
genstudy = False
#Boolean value that determines if the txt with the information of the Gen-particles should be generated or not. Only used when genstudy == True.
gentxt = False


#------------------ HISTOGRAMS ---------------------
if (makeplots == True):

	c1 = ROOT.TCanvas("c1", "Titulo")    # ROOT canvas

	#Creation of empty TH1F objects (empty ROOT histograms)
	plots = histos.histos()   # plots is a list of TH1F objects

	if (checktaus == True):

		plot_tausmalos = TH1F("tausmalos", "tausmalos", 2, 0.0, 3.0)

		plot_PT_tausmenores = TH1F("PT_tausmenores", "PT_tausmenores", 100, 0.0, 1000.0)
		plot_PT_tausmayores = TH1F("PT_tausmayores", "PT_tausmayores", 100, 0.0, 1000.0)
		plot_ETA_tausmenores = TH1F("ETA_tausmenores", "ETA_tausmenores", 100, -5, 5)
		plot_ETA_tausmayores = TH1F("ETA_tausmayores", "ETA_tausmayores", 100, -5, 5)
		plot_PHI_tausmenores = TH1F("PHI_tausmenores", "PHI_tausmenores", 100, -4, 4)
		plot_PHI_tausmayores = TH1F("PHI_tausmayores", "PHI_tausmayores", 100, -4, 4)

		plot_DeltaR_tausmenores = TH1F("DeltaR_tausmenores", "DeltaR_tausmenores", 100, -4, 4)
		plot_DeltaR_tausmayores = TH1F("DeltaR_tausmayores", "DeltaR_tausmayores", 100, -4, 4)
		plot_DeltaPhi_tausmenores = TH1F("DeltaPhi_tausmenores", "DeltaPhi_tausmenores", 100, -4, 4)
		plot_DeltaPhi_tausmayores = TH1F("DeltaPhi_tausmayores", "DeltaPhi_tausmayores", 100, -4, 4)

arrs = []
print("Begin of Data Reading")


#------------- ITERATING THE FILES AND MAKING THE HISTOGRAMS ----------------

for n_signal, signal in enumerate(entries):

	arr1 = []

	#Counter used for saving in txt the information of the Gen-particles of the first 10 events. Only used when gentxt == True
	txt_counter = 0

	#List used to calculate the average total reconstructed PT/PL of the Z' (calculated with tau1 and tau2). Only used when checktaus == True.
	pt_p = []

	#Initializes all the different counters needed for the Gen-particle study in case it must be carried.
	if (genstudy == True):

		#List of PDG-id of first mother of gen pions that have a m1 in all events.
		pion_mother1 = []
		#List of PDG-id of second mother of gen pions that have a m2 in all events.
		pion_mother2 = []
		#List of PDG-id of first daughter of gen taus that have a d1 in all events.
		tau_daughter1 = []
		#List of PDG-id of second daughter of gen taus that have a d2 in all events.
		tau_daughter2 = []
		
		#Counter of times a charged pion's grandmother is a tau/a tau's granddaughter is a charged pion in all events.
		fromtau = 0
		#Number of charged pions with grandmother/tau's with granddaughter in all events.
		withgrand = 0
		#Counters of times the chosen taus from the Z' (tau1 - tau2) match with granddaughter pions from taus (smallest DR) in all events.
		match1 = 0
		match2 = 0

		#Number of good events (match1 and match2 >= 1) in the non-problematic region of tau2: pt(tau2) > 100GeV.
		right_events = 0
		#Number of good events (match1 and match2 >= 1) in the problematic region of tau2: 20GeV < pt(tau2) < 100GeV.
		right_events_p = 0
		#Number of events that pass the cross-cleaning in the non-problematic region of tau2.
		total_events = 0
		#Number of events that pass the cross-cleaning in the problematic region of tau2.
		total_events_p = 0

	if (makeplots == True):
		f = ROOT.TFile(signal + ".root", "recreate")

	for ind in range(1, jobs[n_signal] + 1):

		directory = str("/disco4/SIMULACIONES/Liliana/fully_hadronic/" + signal + "/" + signal + "_" + str(ind) + "/Events/run_01/tag_1_delphes_events.root")
		File = ROOT.TChain("Delphes;1")
		File.Add(directory)
		Number = File.GetEntries()

		print("Signal: " + signal + "_" + str(ind))

		for i in range(Number):
			Entry = File.GetEntry(i)

			#Initializes particles lists.
			jets = []
			bs = []
			METs = []
			taus = []

			if (checktaus == True):
				tausbad = []

			EntryFromBranch_j = File.Jet.GetEntries()
			for j in range(EntryFromBranch_j):

				BTag = File.GetLeaf("Jet.BTag").GetValue(j)
				TauTag = File.GetLeaf("Jet.TauTag").GetValue(j)

				#searches for jets.
				if (BTag != 1 and TauTag != 1):
					jet = TLorentzVector()
					jet_PT, jet_Eta, jet_Phi, jet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
					jet.SetPtEtaPhiM(jet_PT, jet_Eta, jet_Phi, jet_M)
					jets.append(jet)

				#searches for b_jets.
				elif (BTag == 1 and TauTag != 1):
					bjet = TLorentzVector()
					bjet_PT, bjet_Eta, bjet_Phi, bjet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
					bjet.SetPtEtaPhiM(bjet_PT, bjet_Eta, bjet_Phi, bjet_M)
					bs.append(bjet)

				#searches for taus.
				elif (TauTag == 1 and BTag != 1):
					tau = TLorentzVector()
					tau_PT, tau_Eta, tau_Phi, tau_M, tau_Charge = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j),  File.GetLeaf("Jet.Charge").GetValue(j)
					tau.SetPtEtaPhiM(tau_PT, tau_Eta, tau_Phi, tau_M)
					taus.append(tau)

				if (checktaus == True):
					#searches for jets with both tags, bad taus (tauTag and bTag).
					if (TauTag == 1 and BTag == 1):
						tau = TLorentzVector()
						tau_PT, tau_Eta, tau_Phi, tau_M, tau_Charge = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j),  File.GetLeaf("Jet.Charge").GetValue(j)
						tau.SetPtEtaPhiM(tau_PT, tau_Eta, tau_Phi, tau_M)
						tausbad.append(tau)

			# MET (neutrinos).
			Total_MET = 0
			EntryFromBranch_MET = File.MissingET.GetEntries()
			for j in range(EntryFromBranch_MET):
				MET = TLorentzVector()
				MET_PT, MET_Eta, MET_Phi, MET_M  = File.GetLeaf("MissingET.MET").GetValue(j), File.GetLeaf("MissingET.Eta").GetValue(j), File.GetLeaf("MissingET.Phi").GetValue(j), 0.0
				MET.SetPtEtaPhiM(MET_PT, MET_Eta, MET_Phi, MET_M)
				METs.append(MET)
				Total_MET += MET_PT
      

			if (checktaus == True) and (makeplots == True):

				#print(len(taus), len(tausbad))	

				plot_tausmalos.AddBinContent(1, len(taus))
				
				plot_tausmalos.AddBinContent(2, len(tausbad))

				a = plot_tausmalos.GetXaxis()
				a.SetBinLabel(1,"Real taus")
				a.SetBinLabel(2,"Fake taus")


      			#Checks if the event has the minimum theoretical expected number of particles (the count might be higher due to particle-detector effects).
			if (len(jets) >= 4 and len(bs) >= 2 and len(taus) >= 2):

				jets.sort(reverse = True, key=PT)     
				bs.sort(reverse = True, key=PT)
				taus.sort(reverse = True, key=PT)

				#Defines the real jets and bs from the event.
				j1, j2, j3, j4, b1, b2 = jet_reconstruction(jets, bs)
				realjets = [j1, j2, j3, j4]
				realbs = [b1, b2]

				#Makes the cross cleaning, only saving taus with dR not less than 0.3 to any other jet.
				realtaus = cross_cleaning(realjets, realbs, taus)

				#Checks if there are at least two real taus in the event (taus that passed the cross cleaning).
				if (len(realtaus) >= 2):

					tau1, tau2 = tau_reconstruction(realtaus)

					if (checktaus == True):
						pt_p.append(frac_twobody_pTvspL(tau1, tau2))

					#Does all the gen particle studies (matching algorithm, gen txt creation, etc) only if asked to.
					if (genstudy == True):

						genjets = []
						genjetsbad = []
						gen_taus = []
						gen_pions = []
						#Counters of times the chosen taus (tau1 - tau2) match with granddaughter pions (smallest DR) in the event (if it is in the non-problematic region).
						temp_match1 = 0
						temp_match2 = 0

						#Counters of times the chosen taus (tau1 - tau2) match with granddaughter pions (smallest DR) in the event (if it is in the problematic region).
						temp_match1_p = 0
						temp_match2_p = 0
						
						txt_counter += 1
						txt = []

						#Generation jets (to see origin of problematic taus: 20 < PT < 100)
						EntryFromBranch_gen = File.Particle.GetEntries()
						for j in range(EntryFromBranch_gen):

							gen_id = File.GetLeaf("Particle.PID").GetValue(j)


							#~CODE FOR TXT WRITTING.
							#Checks if the txt needs to be written and writes it.
							if (txt_counter < 11) and (txtgen == True):
								str_id = str(gen_id)
								str_PT = str(round(File.GetLeaf("Particle.PT").GetValue(j),1))
								str_ETA = str(round(File.GetLeaf("Particle.Eta").GetValue(j),1))
								str_PHI = str(round(File.GetLeaf("Particle.Phi").GetValue(j),1))
			
								if (File.GetLeaf("Particle.PT").GetValue(j) == 0.0):
									str_PT = "0.000"

								if (File.GetLeaf("Particle.Eta").GetValue(j) == 0.0):
									str_ETA = "0.000"

								m1 = int(File.GetLeaf("Particle.M1").GetValue(j))
								m2 = int(File.GetLeaf("Particle.M2").GetValue(j))

								d1 = int(File.GetLeaf("Particle.D1").GetValue(j))
								d2 = int(File.GetLeaf("Particle.D2").GetValue(j))

								#Checks if the i-th particle has mothers and daughters.
								if (m1 >= 0.0):
									str_m1 = str(File.GetLeaf("Particle.PID").GetValue(m1))
								else: 
									str_m1 = "-"
								
								if (m2 >= 0.0):
									str_m2 = str(File.GetLeaf("Particle.PID").GetValue(m2))
								else: 
									str_m2 = "-"
								

								if (d1 >= 0.0):
									str_d1 = str(File.GetLeaf("Particle.PID").GetValue(d1))
								else: 
									str_d1 = "-"
								
								if (d2 >= 0.0):
									str_d2 = str(File.GetLeaf("Particle.PID").GetValue(d2))
								else: 
									str_d2 = "-"
								
								L = [str_id, str_PT, str_ETA, str_PHI, str_m1, str_m2, str_d1, str_d2]
								txt.append(L)					
							

							#~CODE FOR GEN-PIONS.
							#Checks if the GenParticle is a charged pion or kaon (as around 84% of times a tau decays into a charged pion).
							if ((abs(gen_id) == 211) or (abs(gen_id) == 321)):
							
								genjet = TLorentzVector()
								genjet_PT, genjet_Eta, genjet_Phi, genjet_M  = File.GetLeaf("Particle.PT").GetValue(j), File.GetLeaf("Particle.Eta").GetValue(j), File.GetLeaf("Particle.Phi").GetValue(j), File.GetLeaf("Particle.Mass").GetValue(j)
								genjet.SetPtEtaPhiM(genjet_PT, genjet_Eta, genjet_Phi, genjet_M)

								#Checks if the charged pion is in the problematic region of taus.
								#if ((genjet.Pt() > 15) and (genjet.Pt() < 100)):
									#if (abs(genjet.Eta()) < 2.3):

								genjets.append(genjet)

								#Mother and grandmother index in GenParticle array, < 0 if there is no such particle.
								parton1 = int(File.GetLeaf("Particle.M1").GetValue(j))
								parton2 = int(File.GetLeaf("Particle.M2").GetValue(j))

								if ((parton1 >= 0.0) and (parton2 >= 0.0)):

									#Mother and grandmother ID's of the charged pion.
									parton1_id = File.GetLeaf("Particle.PID").GetValue(parton1)
									parton2_id = File.GetLeaf("Particle.PID").GetValue(parton2)

									pion_mother1.append(abs(parton1_id))
									pion_mother2.append(abs(parton2_id))


							#~CODE FOR GEN-TAUS AND MATCHING ALGORITHM.
							#Checks if the GenParticle is a tau (as around 84% of times a tau decays into a charged pion).
							if (abs(gen_id) == 15):
							
								genjet = TLorentzVector()
								genjet_PT, genjet_Eta, genjet_Phi, genjet_M  = File.GetLeaf("Particle.PT").GetValue(j), File.GetLeaf("Particle.Eta").GetValue(j), File.GetLeaf("Particle.Phi").GetValue(j), File.GetLeaf("Particle.Mass").GetValue(j)
								genjet.SetPtEtaPhiM(genjet_PT, genjet_Eta, genjet_Phi, genjet_M)

								#Checks if the tau is in the problematic region of taus.
								#if ((genjet.Pt() > 15) and (genjet.Pt() < 100)):
									#if (abs(genjet.Eta()) < 2.3):

								gen_taus.append(genjet)
								
								#Daughter and granddaughter index in GenParticle array, < 0 if there is no such particle.
								daughter1 = File.GetLeaf("Particle.D1").GetValue(j)
								daughter2 = File.GetLeaf("Particle.D2").GetValue(j)
								
								if ((daughter1 >= 0.0) and (daughter2 >= 0.0)):

									withgrand += 1
									#Daughter and granddaughter ID's of the gen-tau.
									daughter1_id =  File.GetLeaf("Particle.PID").GetValue(int(daughter1))
									daughter2_id =  File.GetLeaf("Particle.PID").GetValue(int(daughter2))

									tau_daughter1.append(abs(daughter1_id))
									tau_daughter2.append(abs(daughter2_id))

									#Checks if the tau has a charged pion or kaon as a granddaughter.
									if ((abs(daughter2_id) == 211) or (abs(daughter2_id) == 321)):

										fromtau += 1					
										pion = TLorentzVector()
										pion_PT, pion_Eta, pion_Phi, pion_M  = File.GetLeaf("Particle.PT").GetValue(int(daughter2)), File.GetLeaf("Particle.Eta").GetValue(int(daughter2)), File.GetLeaf("Particle.Phi").GetValue(int(daughter2)), File.GetLeaf("Particle.Mass").GetValue(int(daughter2))
										pion.SetPtEtaPhiM(pion_PT, pion_Eta, pion_Phi, pion_M)
										
										gen_pions.append(pion)
							
										#Index and dR of the detector-tau that is closest to the gen-pion.
										best_tau_index = -1
										smallest_DR= 999999.9

										for i in range (len(realtaus)):
											DR = pion.DeltaR(realtaus[i])
											if (DR < smallest_DR):
												smallest_DR = DR
												best_tau_index = i

										
										#Checks if the event is in the non-problematic region.
										if (tau2.Pt() >= 100):
											#Checks if the pion matched with tau1.
											if (best_tau_index == index_tau1):
												match1 += 1
												temp_match1 += 1

											#Checks if the pion matched with tau2.
											elif (best_tau_index == index_tau2):
												match2 += 1
												temp_match2 += 1

										#Checks if the event is in the problematic region.
										else:
											#Checks if the pion matched with tau1.
											if (best_tau_index == index_tau1):
												temp_match1_p += 1
											#Checks if the pion matched with tau2.
											elif (best_tau_index == index_tau2):
												temp_match2_p += 1
									

						#Opens, writes and closes the txt files if necessary.
						if (txt_counter < 11) and (txtgen == True):
							with open(str(signal + "_" + str(txt_counter) + ".txt"), 'w') as g:
								g.write("PDGId      PT           ETA        PHI       Mother   Grandmother    Daughter  Granddaughter")
								g.write('\n')
								for line in txt:
									for i in line:
										g.write(i)
										g.write("        ")
									g.write('\n') 

							g.close()
					
					
						#Checks if the event had the minimun required gen-pions coming from a gen-tau as to count as possible correctly tagged event.
						if (len(gen_pions) >= 2):
						
							#Checks if both taus from the Z' matched with a gen-pion in the non-problematic region.
							if ((temp_match1 >= 1) and (temp_match2 >= 1)):
								right_events += 1

							#Checks if both taus from the Z' matched with a gen-pion in the problematic region.
							if ((temp_match1_p >= 1) and (temp_match2_p >= 1)):
								right_events_p += 1
								
							#Checks if the event was in the non-problematic region.
							if (tau2.Pt() >= 100):
								total_events += 1

							#Checks if the event was in the problematic region.
							if (tau2.Pt() < 100):
								total_events_p += 1
					
	
					#~CONTINUING WITH THE EVENT.				
					#Pt requirement of the different jets (20 for taus is LHC minimum).
					if ((j1.Pt() > 30) and (j2.Pt() > 30) and (j3.Pt() > 30) and (j4.Pt() > 30) and (b1.Pt() > 30) and (j2.Pt() > 30) and (tau1.Pt() > 20) and (tau2.Pt() > 20)):

						#LHC required eta for taus.				
						if ((abs(tau1.Eta()) < 2.3) and (abs(tau2.Eta()) < 2.3)):

							#Minimum possible separation between both taus.
							if (tau1.DeltaR(tau2) > 0.3):

								row = np.array([tau1.Pt(), tau2.Pt(), tau1.DelptaPhi(tau2), tau1.DeltaR(tau2), (tau1 + tau2).M() + Total_MET, (j3 + j4).M(), j1.Pt(), j1.DeltaPhi(j2), j1.DeltaPhi(j3), j1.DeltaPhi(j4)])
            							arr1.append(row)

								#Checks if the plots should be filled and fills them.
								if (makeplots == True):

									#Variables to be plotted.
									variables = [b1.Pt(), b1.Eta(), b1.Phi(), b2.Pt(), b2.Eta(), b2.Phi(), tau1.Pt(), tau1.Eta(), tau1.Phi(), tau2.Pt(), tau2.Eta(), tau2.Phi(), j1.Pt(), j1.Eta(), j1.Phi(), j2.Pt(), j2.Eta(), j2.Phi(), j3.Pt(), j3.Eta(), j3.Phi(), j4.Pt(), j4.Eta(), j4.Phi(), tau1.DeltaR(tau2), b1.DeltaR(b2), j1.DeltaR(j2), j1.DeltaR(j3), j1.DeltaR(j4), j2.DeltaR(j3), j2.DeltaR(j4), j3.DeltaR(j4), tau1.DeltaPhi(tau2), b1.DeltaPhi(b2), j1.DeltaPhi(j2), j1.DeltaPhi(j3), j1.DeltaPhi(j4), j2.DeltaPhi(j3), j2.DeltaPhi(j4), j3.DeltaPhi(j4) , abs(tau1.Pt() - tau2.Pt()), abs(b1.Pt() - b2.Pt()), abs(j1.Pt() - j2.Pt()), abs(j1.Pt() - j3.Pt()), abs(j1.Pt() - j4.Pt()), abs(j2.Pt() - j3.Pt()), abs(j2.Pt() - j4.Pt()), abs(j3.Pt() - j4.Pt()) , (tau1 + tau2).M() + Total_MET, (j1+j2).M(), (j1+j3).M(), (j1+j4).M(), (j2+j3).M(), (j2+j4).M(), (j3+j4).M(), (j1+j2+b1).M(), (j1+j3+b1).M(), (j1+j4+b1).M(), (j2+j3+b1).M(), (j2+j4+b1).M(), (j3+j4+b1).M(), (j1+j2+b2).M(), (j1+j3+b2).M(), (j1+j4+b2).M(), (j2+j3+b2).M(), (j2+j4+b2).M(), (j3+j4+b2).M(), j1.Pt() + j2.Pt() + j3.Pt() + j4.Pt() + b1.Pt() + b2.Pt() + tau1.Pt() + tau2.Pt() + Total_MET, taus[0].Pt(), taus[1].Pt(), frac_twobody_pTvspL(tau1, tau2), tau1.P(), tau2.P()]

									
									for i in range(len(plots)):
								  		histos_fill(plots[i], variables[i])
									
									if (checktaus == True):
										if (tau2.Pt() < 80.0):
											histos_fill(plot_PT_tausmenores, tau2.Pt())
											histos_fill(plot_ETA_tausmenores, tau2.Eta())
											histos_fill(plot_PHI_tausmenores, tau2.Phi())

											histos_fill(plot_DeltaR_tausmenores, tau1.DeltaR(tau2))
											histos_fill(plot_DeltaPhi_tausmenores, tau1.DeltaPhi(tau2))
										
										if (tau2.Pt() > 80.0):
											histos_fill(plot_PT_tausmayores, tau2.Pt())
											histos_fill(plot_ETA_tausmayores, tau2.Eta())
											histos_fill(plot_PHI_tausmayores, tau2.Phi())

											histos_fill(plot_DeltaR_tausmayores, tau1.DeltaR(tau2))
											histos_fill(plot_DeltaPhi_tausmayores, tau1.DeltaPhi(tau2))
	
	
	arrs.append(arr1)

	if (genstudy == True):		
							
		c = Counter(pion_mother1)
		c.most_common(1)
		print ("pions most common mother: ", len(pion_mother1), c.most_common(1))

		c2 = Counter(pion_mother2)
		c2.most_common(1)
		print ("pions most common grandmother: ", len(pion_mother2), c2.most_common(1))

		c3 = Counter(tau_daughter1)
		c3.most_common(1)
		print ("taus most common daughter: ", len(tau_daughter1), c3.most_common(1))
		
		c4 = Counter(tau_daughter2)
		c4.most_common(1)
		print ("taus most common granddaughter: ", len(tau_daughter2), c4.most_common(1))
		
		print ("taus with granddaughter: ", withgrand,"taus with pion as a granddaughter: ", fromtau, "taus1 that match with a pion: ", match1, "taus2 that match with a pion: ", match2, "total events that pass cross cleaning in the non-problematic region: ", total_events, "events with tau1 and tau2 matching in the non-problematic region: ", right_events, "total events that pass cross cleaning in the problematic region: ", total_events_p, "events with tau1 and tau2 matching in the problematic region: ", right_events_p)

	if (checktaus == True):
		#Average total reconstructed PT/PL of the Z' (calculated with tau1 and tau2).
		print (sum(pt_p)/len(pt_p))


	if (makeplots == True):

		#Drawing the histograms.
		for plot in plots:
			histos_Draw(plot)

		if (checktaus == True):
			plot_tausmalos.Draw('HISTOS')

			plot_PT_tausmenores.Draw('HISTOS')
			plot_PT_tausmayores.Draw('HISTOS')
			plot_ETA_tausmenores.Draw('HISTOS')
			plot_ETA_tausmayores.Draw('HISTOS')
			plot_PHI_tausmenores.Draw('HISTOS')
			plot_PHI_tausmayores.Draw('HISTOS')

			plot_DeltaR_tausmenores.Draw('HISTOS')
			plot_DeltaR_tausmayores.Draw('HISTOS')
			plot_DeltaPhi_tausmenores.Draw('HISTOS')
			plot_DeltaPhi_tausmayores.Draw('HISTOS')


		# Updating the canvas.
		c1.Update()


		# Writing the histograms.
		for plot in plots:
			histos_Write(plot)

		if (checktaus == True):
			plot_tausmalos.Write()

			plot_PT_tausmenores.Write()
			plot_PT_tausmayores.Write()
			plot_ETA_tausmenores.Write()
			plot_ETA_tausmayores.Write()
			plot_PHI_tausmenores.Write()
			plot_PHI_tausmayores.Write()

			plot_DeltaR_tausmenores.Write()
			plot_DeltaR_tausmayores.Write()
			plot_DeltaPhi_tausmenores.Write()
			plot_DeltaPhi_tausmayores.Write()

		# Closing the ROOT file where the histos were saved.
		f.Close()

		  
		# Reseting the TH1F objects for its use in the next signal or bkg file.
		for plot in plots:
			histos_Reset(plot)

		if (checktaus == True):
			plot_tausmalos.Reset('ICESM')

			plot_PT_tausmenores.Reset('ICESM')
			plot_PT_tausmayores.Reset('ICESM')
			plot_ETA_tausmenores.Reset('ICESM')
			plot_ETA_tausmayores.Reset('ICESM')
			plot_PHI_tausmenores.Reset('ICESM')
			plot_PHI_tausmayores.Reset('ICESM')

			plot_DeltaR_tausmenores.Reset('ICESM')
			plot_DeltaR_tausmayores.Reset('ICESM')
			plot_DeltaPhi_tausmenores.Reset('ICESM')
			plot_DeltaPhi_tausmayores.Reset('ICESM')


print("End of data reading")

print("Begin of Data Preparation")


arrs = np.array(arrs)
bkg1 = np.array(arrs[0])
bkg2 = np.array(arrs[1])
bkg3 = np.array(arrs[2])
signal1, signal1_vsize = np.array(arrs[3]), np.shape(np.array(arrs[3]))[0]
signal2, signal2_vsize = np.array(arrs[4]), np.shape(np.array(arrs[4]))[0]
signal3, signal3_vsize = np.array(arrs[5]), np.shape(np.array(arrs[5]))[0]

pred1 = np.concatenate((signal1, bkg1[:int(signal1_vsize/3.0), :], bkg2[:int(signal1_vsize/3.0), :], bkg3[:signal1_vsize - 2*int(signal1_vsize/3.0), :]), axis=0)
pred2 = np.concatenate((signal2, bkg1[:int(signal2_vsize/3.0), :], bkg2[:int(signal2_vsize/3.0), :], bkg3[:signal2_vsize - 2*int(signal2_vsize/3.0), :]), axis=0)
pred3 = np.concatenate((signal3, bkg1[:int(signal3_vsize/3.0), :], bkg2[:int(signal3_vsize/3.0), :], bkg3[:signal3_vsize - 2*int(signal3_vsize/3.0), :]), axis=0)


labels1 = np.zeros(np.shape(pred1)[0])
labels2 = np.zeros(np.shape(pred2)[0])
labels3 = np.zeros(np.shape(pred3)[0])

labels1[:signal1_vsize] = 1
labels2[:signal2_vsize] = 1
labels3[:signal3_vsize] = 1


print("End of Data Preparation")

print("Executing train_test_split...")
trainPredictors_m350, testPredictors_m350, trainLabels_m350, testLabels_m350 = train_test_split(pred1, labels1, test_size=0.25)
trainPredictors_m1000, testPredictors_m1000, trainLabels_m1000, testLabels_m1000 = train_test_split(pred2, labels2, test_size=0.25)
trainPredictors_m3000, testPredictors_m3000, trainLabels_m3000, testLabels_m3000 = train_test_split(pred3, labels3, test_size=0.25)

testPredictors = [testPredictors_m350, testPredictors_m1000, testPredictors_m3000]
testLabels = [testLabels_m350, testLabels_m1000, testLabels_m3000]


logreg_m350_model = bcml_model(make_pipeline(StandardScaler(), LogisticRegression()))
logreg_m350_model.fit(trainPredictors_m350, trainLabels_m350)

logreg_m1000_model = bcml_model(make_pipeline(StandardScaler(), LogisticRegression()))
logreg_m1000_model.fit(trainPredictors_m1000, trainLabels_m1000)

logreg_m3000_model = bcml_model(make_pipeline(StandardScaler(), LogisticRegression()))
logreg_m3000_model.fit(trainPredictors_m3000, trainLabels_m3000)


masses, sig_css, bg_css = get_elijah_ttbarzp_cs()

masses = [350, 1000, 3000]
sig_css = [0.001813, 0.0002066, 3.114E-6]

conv = 10**15/10**12 # conv * lumi (in fb-1)*crossx (im pb) = # of events
lumi = 3000

signal_yields = [conv*lumi*sig_cs for sig_cs in sig_css]
background_yields = [conv*lumi*bg_cs for bg_cs in bg_css]
background_yield = sum(background_yields)

bg_labels = [r"No Higgs B.g. ($pp \; \to \; t\overline{t}b\overline{b} \; \backslash \; h \; \to \; bjj \; \overline{b}\ell\nu \;  b\overline{b}$)",
             r"Higgs B.g. ($pp \; \to \; t\overline{t}h \; \to \; bjj \; \overline{b}\ell\nu \;  b\overline{b}$)",
             r"4 Tops B.g. ($pp \; \to \; t\overline{t}t\overline{t} \; \to \; \overline{b}\ell\nu \;  b\overline{b} \; ...$)"]


sample_masses = np.linspace(masses[0], masses[-1], 500)
sigs_f = scipy.interpolate.interp1d(masses, np.log10(sig_css), kind='quadratic')
new_settings = get_settings()
new_settings['figure.figsize'] = (15,10)
with plt.rc_context(new_settings):
    plt.plot(sample_masses, [10**sigs_f(m) for m in sample_masses], 
             label=r"Signal ($pp \; \to \; t\overline{t}Z' \; \to \; bjj \; \overline{b}\ell\nu \;  b\overline{b}$)")
    for bg_cs, label in zip([bg_css[2], bg_css[0], bg_css[1]], bg_labels):
        plt.plot([350, 3000], [bg_cs]*2, label=r"{}".format(label))
    plt.ylabel('Cross Section (pb)');
    plt.xlabel(r"$Z'$ mass $m_{Z'}$ (GeV)");
    plt.xlim(350, 3000)
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(0,0), loc="lower left")
    plt.savefig('crossx.png', dpi=200, bbox_inches='tight')
    plt.clf()


filename_bgs = entries[:3]

for B, name in zip(background_yields, filename_bgs):
    print(f"Background {name} constitutes {round(B/background_yield,3)}% of the total background cross section.")
print(f"\nTotal backgroundcross section: {sum(bg_css)} pb")


def getTime(form='s'):
    seconds = datetime.now().day * 86400 + datetime.now().hour * 3600 + datetime.now().minute * 60 + datetime.now().second
    if form == 's':
        return seconds
    elif form == 'm':
        return round(seconds/60,3)
    elif form == 'h':
        return round(seconds/3600,3)

logreg_models = [logreg_m350_model, logreg_m1000_model, logreg_m3000_model]


time_before = getTime()

logreg_sigs_nosepbg = [model.significance(signal_yield, background_yield, tpr=model.tpr(testLabels[i], preds=testPredictors[i]), fpr=model.tpr(testLabels[i], preds=testPredictors[i]), sepbg=False) for i, (model, signal_yield) in enumerate(zip(logreg_models, signal_yields))]

for sig, mass in zip(logreg_sigs_nosepbg, masses):
    print(f"Z' mass = {mass} GeV --> logistic regression gives significance of {round(sig, 3)} sigma", "w/o separated backgrounds")
time_after = getTime()
print(f"(Runtime: {time_after - time_before} seconds)")

"""
zp_cs = cross_section_helper(masses, sig_css, bg_css, mass_units='GeV')

max_mass = zp_cs.absolute_max_mass_sens()
logreg_sigs_nosepbg_f = scipy.interpolate.interp1d(masses, np.log10(logreg_sigs_nosepbg), kind='quadratic')

with plt.rc_context(get_settings()):
    plt.plot([max_mass for i in range(2)], np.logspace(-2, 1.5, 2), label=r'Theoretical $5\sigma$ Limit', c='black')
    plt.plot(sample_masses, [5 for m in sample_masses], label=r'5$\sigma$', c='blue')
    plt.plot(sample_masses, [1.645 for m in sample_masses], label=r'90% confidence', c='cornflowerblue')
    plt.plot(
        sample_masses, [10**logreg_sigs_nosepbg_f(m) for m in sample_masses], label='Logistic Regression (No b.g. sep.)', 
        c='red')
    plt.ylabel('Signal Significance');
    plt.xlabel(r"$Z'$ mass $m_{Z'}$ (GeV)");
    plt.yscale('log')
    plt.legend()
"""

background_names = entries[:3]
signal_names_2 = ["m(Z')=350 GeV", "m(Z')=1000 GeV", "m(Z')=3000 GeV"]

with plt.rc_context(get_settings()):
    for i, model in enumerate(logreg_models):
        plt.figure(i)
        bin_edges, sig_bins, bg_bins = model.predict_hist(testPredictors[i], testLabels[i], num_bins=40, sepbg=False)
        plt.bar(bin_edges[:-1], sig_bins, width=np.diff(bin_edges), alpha=0.75, align="edge", label="Signal")
        plt.bar(bin_edges[:-1], bg_bins, width=np.diff(bin_edges), alpha=0.75, align="edge", label="Background")
        plt.yscale('log')
        plt.title('Threshold Optimization for {}'.format(signal_names_2[i]))
        plt.gca().set_ylim(10**-3, 10**2)
        plt.gca().set_xlim(0, 1)
        plt.legend()
        plt.savefig(f'log_reg_m{masses[i]}GeV_hist.png', dpi=200, bbox_inches='tight')
        plt.clf()

time_before = getTime(form='m')
logreg_opt_results_nosepbg = [model.best_threshold(signal_yield, background_yield, testPredictors[i], testLabels[i], sepbg=False) for i, (model, signal_yield) in enumerate(zip(logreg_models, signal_yields))]

logreg_opt_sigs_nosepbg = [result[1] for result in logreg_opt_results_nosepbg]
for sig, mass in zip(logreg_opt_sigs_nosepbg, masses):
    print(f"Z' mass = {mass} GeV --> logistic regression gives optimized sign. of {round(sig, 3)} sigma",
          "w/o separated backgrounds")
time_after = getTime(form='m')
print(f"(Runtime: {round(time_after - time_before,3)} minutes)")


with plt.rc_context(get_settings()):
    for i, (model, result) in enumerate(zip(logreg_models, logreg_opt_results_nosepbg)):
        plt.figure(i)
        bin_edges, sig_bins, bg_bins = model.predict_hist(testPredictors[i], testLabels[i], num_bins=40, sepbg=False)
        plt.bar(bin_edges[:-1], sig_bins, width=np.diff(bin_edges), alpha=0.75, align="edge", label="Signal")
        plt.bar(bin_edges[:-1], bg_bins, width=np.diff(bin_edges), alpha=0.75, align="edge", label="Background")
        plt.plot([0.5]*2, [10**-3, 10**2], label='Default Threshold', c='tab:red')
        plt.plot([result[0]]*2, [10**-3, 10**2], label='Optimized Threshold', c='tab:green')
        plt.yscale('log')
        plt.title('Threshold Optimization for {}'.format(signal_names_2[i]))
        plt.gca().set_ylim(10**-3, 10**2)
        plt.gca().set_xlim(0, 1)
        plt.legend()
        plt.savefig(f'log_reg_m{masses[i]}GeV_opt_hist', dpi=200, bbox_inches='tight')

##########################

time_before = getTime(form='m')
xgradboost_m350G_model = bcml_model(
    make_pipeline(StandardScaler(), 
                  XGBClassifier(n_estimators=250, max_depth=7, learning_rate=0.1, nthread=4, use_label_encoder=False)));
xgradboost_m350G_model.fit(trainPredictors_m350, trainLabels_m350)
time_after = getTime(form='m')
print(f"(Runtime: {round(time_after - time_before,3)} minutes)")

time_before = getTime(form='m')
xgradboost_m1000G_model = bcml_model(
    make_pipeline(StandardScaler(), 
                  XGBClassifier(n_estimators=250, max_depth=7, learning_rate=0.1, nthread=4, use_label_encoder=False)));
xgradboost_m1000G_model.fit(trainPredictors_m1000, trainLabels_m1000)
time_after = getTime(form='m')
print(f"(Runtime: {round(time_after - time_before,3)} minutes)")

time_before = getTime(form='m')
xgradboost_m3000G_model = bcml_model(
    make_pipeline(StandardScaler(), 
                  XGBClassifier(n_estimators=250, max_depth=7, learning_rate=0.1, nthread=4, use_label_encoder=False)));
xgradboost_m3000G_model.fit(trainPredictors_m3000, trainLabels_m3000)
time_after = getTime(form='m')
print(f"(Runtime: {round(time_after - time_before,3)} minutes)")

xgradboost_models = [
    xgradboost_m350G_model, xgradboost_m1000G_model, xgradboost_m3000G_model]

time_before = getTime()
xgradboost_sigs_nosepbg = [model.significance(signal_yield, background_yield, tpr=model.tpr(testLabels[i], preds=testPredictors[i]), fpr=model.tpr(testLabels[i], preds=testPredictors[i]), sepbg=False) for i, (model, signal_yield) in enumerate(zip(xgradboost_models, signal_yields))]

for sig, mass in zip(xgradboost_sigs_nosepbg, masses):
    print(f"Z' mass = {mass} GeV --> gradient boosting gives significance of {round(sig, 3)} sigma",
          "w/o separated backgrounds")
time_after = getTime()
print(f"(Runtime: {time_after - time_before} seconds)")


time_before = getTime(form='m')
xgradboost_opt_results_nosepbg = [model.best_threshold(signal_yield, background_yield, testPredictors[i], testLabels[i], sepbg=False) for i, (model, signal_yield) in enumerate(zip(xgradboost_models, signal_yields))]

xgradboost_opt_sigs_nosepbg = [result[1] for result in xgradboost_opt_results_nosepbg]
for sig, mass in zip(xgradboost_opt_sigs_nosepbg, masses):
    print(f"Z' mass = {mass} GeV --> gradient boosting forests gives optimized sign. of {round(sig, 3)} sigma",
          "w/o separated backgrounds")
time_after = getTime(form='m')
print(f"(Runtime: {round(time_after - time_before,3)} minutes)")

for mass, result, sig_cs in zip(masses, xgradboost_opt_results_nosepbg, sig_css):
    print(f"for mass = {mass} GeV, thresh = {round(result[0],6)}, tpr = {round(result[2],5)} and fpr = {round(result[3],9)}",
          f"so yield is {int(lumi * sig_cs * conv * result[2])} and sig. sign. is {round(result[1],3)}")

xgradboost_sigs_nosepbg_f = scipy.interpolate.interp1d(masses, np.log10(xgradboost_sigs_nosepbg), kind='quadratic')
xgradboost_opt_sigs_nosepbg_f = scipy.interpolate.interp1d(masses, np.log10(xgradboost_opt_sigs_nosepbg), kind='quadratic')

masses, sig_css, bg_css = get_elijah_ttbarzp_cs()
zp_cs = cross_section_helper(masses, sig_css, bg_css, mass_units='GeV')
max_mass = zp_cs.absolute_max_mass_sens()

with plt.rc_context(get_settings()):
    plt.plot([max_mass for i in range(2)], np.logspace(-2, 1.5, 2), label=r'Theoretical $5\sigma$ Limit', c='black')
    plt.plot(sample_masses, [5 for m in sample_masses], label=r'5$\sigma$', c='blue')
    plt.plot(sample_masses, [1.645 for m in sample_masses], label=r'90% confidence', c='cornflowerblue')
    plt.plot(
        sample_masses, [10**xgradboost_sigs_nosepbg_f(m) for m in sample_masses], 
        label='Gradient Boosting (No b.g. sep.)', 
        c='tab:purple')
    plt.plot(
        sample_masses, [10**xgradboost_opt_sigs_nosepbg_f(m) for m in sample_masses], 
        label='Gradient Boosting (No b.g. sep., opt. thresh.)', 
        c='tab:purple', linestyle='--')
    plt.ylabel('Signal Significance');
    plt.xlabel(r"$Z'$ mass $m_{Z'}$ (GeV)");
    plt.yscale('log')
    plt.legend()
    plt.savefig("grad_boost_sigs.png", dpi=200, bbox_inches='tight')

logreg_opt_sigs_nosepbg_f = scipy.interpolate.interp1d(masses, np.log10([float(item) for item in logreg_opt_results_nosepbg]), kind='quadratic')


with plt.rc_context(get_settings()):
    plt.plot([max_mass for i in range(2)], np.logspace(-2.5, 2, 2), label=r'Theoretical $5\sigma$ Limit', c='black')
    plt.plot(sample_masses, [5 for m in sample_masses], label=r'5$\sigma$', c='blue')
    plt.plot(sample_masses, [1.645 for m in sample_masses], label=r'90% confidence', c='cornflowerblue')
    plt.plot(
        sample_masses, [10**logreg_opt_sigs_nosepbg_f(m) for m in sample_masses], 
        label='Logistic Regression', 
        c='red', linestyle='-.')
    plt.plot(
        sample_masses, [10**xgradboost_opt_sigs_nosepbg_f(m) for m in sample_masses], 
        label='Gradient Boosting', 
        c='tab:purple', linestyle='-.')
    plt.ylabel('Signal Significance');
    plt.xlabel(r"$Z'$ mass $m_{Z'}$ (GeV)");
    plt.yscale('log')
    plt.legend()
    plt.savefig("all_model_sigs.png", dpi=200, bbox_inches='tight')

















